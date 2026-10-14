#!/usr/bin/env python3
"""Rewrite CCEX inputs written with the input keys before 2026-10-13 (A47) in the new keys.

    python3 convert_input_keys.py FILE...            rewrite each file in place and list what changed
    python3 convert_input_keys.py --check FILE...    list only; exit 1 when a file still has old keys
    python3 convert_input_keys.py --stdout FILE      print the new text, leave the file alone

The table of old and new names is docs/user/InputFile.md; (internal notes) has the decisions.
Comments (# or ! to the end of a line, outside strings) stay: a key that moves into a block takes the comment lines
above it and the comment after it on its line. A value is copied as written unless its form changes (method, a 0/1
switch, hf_readmode, ...). An object nothing changed in is left byte for byte. A file in the new form is left alone.
A key written twice (an old key beside its new name, or two spellings that differ in case only) is kept once when the
values are the same; different values stop the file before anything is written. Tests: tests/ (python3 -m unittest).
"""
import json
import os
import shutil
import sys
import tempfile

WS = ' \t\r\n'


class Node:
    def __init__(self, text, start, end):
        self.text, self.start, self.end, self.dirty = text, start, end, False

    def raw(self):
        return self.text[self.start:self.end]


class Scalar(Node):
    def __init__(self, raw_text):
        super().__init__(raw_text, 0, len(raw_text))

    def render(self):
        return self.text


class Member:
    """pre: blank lines and comments before the key; tail: what follows the value on its line (a comment)."""
    def __init__(self, pre, keyraw, mid, value, tail='', shift=0):
        self.pre, self.keyraw, self.mid, self.value, self.tail, self.shift = pre, keyraw, mid, value, tail, shift

    @property
    def key(self):
        return json.loads(self.keyraw)

    def indent(self):
        return self.pre.rsplit('\n', 1)[-1] if '\n' in self.pre else None


class Container(Node):
    open_ch = close_ch = ''

    def __init__(self, text, start, end, items, close):
        super().__init__(text, start, end)
        self.items, self.close = items, close

    def render(self):
        if not self.dirty and not any(changed(m.value) for m in self.items):
            return self.raw()
        out = [self.open_ch]
        for i, m in enumerate(self.items):
            body = render_value(m.value, m.shift)
            last = i == len(self.items) - 1
            out.append(m.pre + (m.keyraw + m.mid if m.keyraw is not None else '') + body + ('' if last else ',') + m.tail)
        out.append(self.close + self.close_ch)
        return ''.join(out)


class Obj(Container):
    open_ch, close_ch = '{', '}'

    def find(self, key):
        for i, m in enumerate(self.items):
            if m.key.lower() == key.lower():
                return i, m
        return None, None

    def get(self, key):
        return self.find(key)[1]

    def take(self, key):
        """Remove the member. The comment lines right above it go with it; what is above them (blank lines, comments
        that belong to the lines before) stays where it was, in front of the member that comes next."""
        i, m = self.find(key)
        if m is not None:
            del self.items[i]
            self.dirty = True
            lines = m.pre.split('\n')
            k = len(lines) - 1
            while k > 0 and lines[k - 1].strip():
                k -= 1
            keep, m.pre = '\n'.join(lines[:k]), '\n' + '\n'.join(lines[k:])
            if keep.strip():
                if i < len(self.items):
                    nxt = self.items[i]
                    nxt.pre = keep + nxt.pre
                    nxt.left = keep + getattr(nxt, 'left', '')
                else:
                    self.close = keep + self.close
            elif keep and i < len(self.items):
                self.items[i].blank = keep          # blank lines that a member put in its place takes over
        return i, m

    def drop(self, key):
        """Remove the member for good; its comment lines stay where they were."""
        i, m = self.find(key)
        if m is None:
            return None
        return self.drop_at(i)

    def drop_at(self, i):
        m = self.items[i]
        del self.items[i]
        self.dirty = True
        body = m.pre.rsplit('\n', 1)[0] if '\n' in m.pre else ''
        if body.strip() or (body and i < len(self.items) and not self.items[i].pre.startswith('\n\n')):
            if i < len(self.items):
                self.items[i].pre = body + self.items[i].pre
            else:
                self.close = body + self.close
        return m

    def insert(self, i, member):
        """A new member at i: the text another removal left in front of the member there goes in front of it."""
        if i < len(self.items):
            nxt = self.items[i]
            left = getattr(nxt, 'left', '')
            if left and nxt.pre.startswith(left):
                nxt.pre, nxt.left = nxt.pre[len(left):], ''
                member.pre = left + member.pre
            elif getattr(nxt, 'blank', '') and '\n' in member.pre:
                member.pre, nxt.blank = nxt.blank + member.pre, ''
        self.items.insert(i, member)
        self.dirty = True


class Arr(Container):
    open_ch, close_ch = '[', ']'


def changed(node):
    return isinstance(node, Container) and (node.dirty or any(changed(m.value) for m in node.items))


def render_value(node, shift=0):
    text = node.render()
    if shift and '\n' in text:
        lines = text.split('\n')
        if shift > 0:
            lines = [lines[0]] + [(' ' * shift + l) if l.strip() else l for l in lines[1:]]
        else:
            lines = [lines[0]] + [l[min(-shift, len(l) - len(l.lstrip(' '))):] for l in lines[1:]]
        text = '\n'.join(lines)
    return text


class Parser:
    def __init__(self, text):
        self.t = text

    def ws(self, i):
        t = self.t
        while i < len(t):
            if t[i] in WS:
                i += 1
            elif t[i] in '#!':
                while i < len(t) and t[i] != '\n':
                    i += 1
            else:
                break
        return i

    def string_end(self, i):
        t = self.t
        assert t[i] == '"', f'expected a string at {i}'
        i += 1
        while t[i] != '"':
            i += 2 if t[i] == '\\' else 1
        return i + 1

    def tail(self, i, close):
        """After a value: (tail text, index where the next member's text starts)."""
        t = self.t
        j = self.ws(i)
        if t[j] == ',':
            k = j + 1                       # a comment after the comma on its line belongs to this member;
            while k < len(t) and t[k] in ' \t':   # the spaces before the next key on the same line, to that key
                k += 1
            if k < len(t) and t[k] in '#!':
                while k < len(t) and t[k] != '\n':
                    k += 1
            else:
                k = j + 1
            return t[i:j] + t[j + 1:k], k, False
        assert t[j] == close, f'expected , or {close} at {j}'
        nl = t.find('\n', i, j)
        if nl < 0:
            return t[i:j], j, True
        return t[i:nl], nl, True

    def value(self, i):
        t = self.t
        c = t[i]
        if c == '{' or c == '[':
            close = '}' if c == '{' else ']'
            items, k = [], i + 1
            while True:
                j = self.ws(k)
                if t[j] == close:
                    node = (Obj if c == '{' else Arr)(t, i, j + 1, items, t[k:j])
                    return node, j + 1
                if c == '{':
                    ke = self.string_end(j)
                    colon = self.ws(ke)
                    assert t[colon] == ':', f'expected : at {colon}'
                    v0 = self.ws(colon + 1)
                    val, e = self.value(v0)
                    tail, nxt, last = self.tail(e, close)
                    items.append(Member(t[k:j], t[j:ke], t[ke:v0], val, tail))
                else:
                    val, e = self.value(j)
                    tail, nxt, last = self.tail(e, close)
                    items.append(Member(t[k:j], None, '', val, tail))
                k = nxt
                if last:
                    j = self.ws(k)
                    node = (Obj if c == '{' else Arr)(t, i, j + 1, items, t[k:j])
                    return node, j + 1
        if c == '"':
            e = self.string_end(i)
        else:
            e = i
            while e < len(t) and t[e] not in WS + ',}]#!':
                e += 1
        node = Node(t, i, e)
        node.render = node.raw
        return node, e


def plain(node):
    """The Python value of a node (comments stripped)."""
    return json.loads(strip_comments(node.render()))


def strip_comments(text):
    out, in_str, esc, i = [], False, False, 0
    while i < len(text):
        c = text[i]
        if in_str:
            out.append(c)
            if esc:
                esc = False
            elif c == '\\':
                esc = True
            elif c == '"':
                in_str = False
        elif c in '#!':
            while i < len(text) and text[i] != '\n':
                i += 1
            continue
        else:
            if c == '"':
                in_str = True
            out.append(c)
        i += 1
    return ''.join(out)


class Converter:
    def __init__(self, text):
        self.text = text
        p = Parser(text)
        s = p.ws(0)
        self.root, e = p.value(s)
        if not isinstance(self.root, Obj):
            raise ValueError('the input is not a JSON object')
        self.prefix, self.suffix = text[:s], text[e:]
        self.log = []
        self.unit = self.detect_unit()
        first = self.root.items[0] if self.root.items else None
        self.mid = first.mid if first is not None else ': '
        self.width = len(first.keyraw) if first is not None else 0      # keys padded to a column, as the file pads them

    def detect_unit(self):
        top = self.member_indent(self.root) or ''
        steps = []                              # how much deeper the entries of a block or a list are written
        for m in self.root.items:
            if isinstance(m.value, Container) and m.value.items:
                inner = m.value.items[0].indent()
                if inner is not None and len(inner) > len(top):
                    steps.append(len(inner) - len(top))
        return ' ' * min(steps) if steps and min(steps) <= 8 else '    '

    @staticmethod
    def member_indent(obj):
        for m in obj.items:
            ind = m.indent()
            if ind is not None:
                return ind
        return None

    def indent_of(self, obj, parent_indent):
        ind = self.member_indent(obj)
        return ind if ind is not None else parent_indent + self.unit

    def note(self, text):
        self.log.append(text)

    # ---- building ---------------------------------------------------------------------------------------------
    @staticmethod
    def aligned(mid, width, keyraw):
        """mid for keyraw, keeping the colon in the column it was in (key + spaces = width)."""
        lead = len(mid) - len(mid.lstrip(' '))
        if lead == 0:
            return mid
        return ' ' * max(1, width + lead - len(keyraw)) + mid.lstrip(' ')

    def new_member(self, key, value, indent, moved=None, from_indent=None):
        """A member for key at indent: a new value (text or node), or the value of a moved member with its comments."""
        if moved is not None:
            comments = [l.strip() for l in moved.pre.split('\n') if l.strip()]
            pre = ''.join('\n' + indent + c for c in comments) + '\n' + indent
            # the moved value keeps its own lines, shifted to the new depth; a new value is built at its depth
            shift = len(indent) - len(from_indent) if (from_indent is not None and value is None) else 0
            node = moved.value if value is None else value
            if isinstance(node, str):
                node = Scalar(node)
            keyraw = json.dumps(key)
            return Member(pre, keyraw, self.aligned(moved.mid, len(moved.keyraw), keyraw), node, moved.tail, shift)
        node = Scalar(value) if isinstance(value, str) else value
        keyraw = json.dumps(key)
        return Member('\n' + indent, keyraw, self.aligned(self.mid, self.width, keyraw), node, '')

    def rename(self, member, key):
        keyraw = json.dumps(key)
        member.mid = self.aligned(member.mid, len(member.keyraw), keyraw)
        member.keyraw = keyraw

    def new_block(self, indent):
        """An empty object whose members will sit at indent + unit and whose } sits at indent."""
        o = Obj('', 0, 0, [], '\n' + indent)
        o.dirty = o.new = True
        return o

    def child_block(self, parent, parent_indent, key, position=None):
        """parent[key] as an object, made (at position, or at the end) when absent."""
        m = parent.get(key)
        if m is not None:
            if not isinstance(m.value, Obj):
                raise ValueError(f'"{key}" is not an object')
            return m.value
        o = self.new_block(parent_indent)
        mem = self.new_member(key, o, parent_indent)
        if position is None:
            parent.items.append(mem)
        else:
            parent.insert(position, mem)
        parent.dirty = True
        return o

    def put(self, block, block_indent, key, moved=None, value=None, from_indent=None):
        ind = self.indent_of(block, block_indent)
        block.items.append(self.new_member(key, value, ind, moved, from_indent))
        block.dirty = True

    # ---- the conversion -----------------------------------------------------------------------------------------
    def convert(self):
        R = self.root
        top = self.member_indent(R) or ''
        self.method_shape(R)
        self.radii(R)
        self.renames(R)
        self.savemode(R, top)
        self.blocks(R, top)
        self.qubit(R, top)
        self.hamiltonian(R, top)
        self.defect(R, top)
        self.single_keys(R, '')
        self.compact(R)
        return self.prefix + R.render() + self.suffix

    def single_keys(self, node, path):
        """One member per key, case ignored as CCEX ignores it (review 20261013-180049 C1). An old key whose new name is
        already written (npulse beside pulse.count) gave the key twice, and CCEX reads the first, Python json the last:
        the same value is kept once, different values refuse the file before anything is written."""
        if isinstance(node, Obj):
            i = 0
            while i < len(node.items):
                m = node.items[i]
                first = next(x for x in node.items if x.key.lower() == m.key.lower())
                if first is m:
                    i += 1
                    continue
                a, b = plain(first.value), plain(m.value)
                if json.dumps(a, sort_keys=True) != json.dumps(b, sort_keys=True):     # true is not 1 here
                    raise ValueError(f'{path}{first.key} is written twice with different values ({json.dumps(a)} and '
                                     f'{json.dumps(b)}): an old key and its new name were both in the input, or a key in '
                                     f'two cases; keep one')
                node.drop_at(i)                 # its comment lines stay; the comment after it goes above the one kept
                if m.tail.strip():
                    ind = first.indent()
                    if ind is not None:
                        first.pre = first.pre[:len(first.pre) - len(ind)] + ind + m.tail.strip() + '\n' + ind
                    else:
                        self.note(f'{path}{m.key}: its comment {m.tail.strip()!r} was not kept (the block is on one line)')
                self.note(f'{path}{m.key}: the same value was written twice, kept once')
        for n, m in enumerate(node.items):
            if isinstance(m.value, Container):
                self.single_keys(m.value, f'{path}{m.key}.' if m.keyraw is not None else f'{path[:-1]}[{n}].')

    def compact(self, obj):
        """A new block of a few plain values on one line: "pulse": {"count": 1}. The comment lines that came with its
        values go above it; a block holding a block, or a value with a comment after it, stays as it is."""
        for m in obj.items:
            v = m.value
            if not isinstance(v, Container):
                continue
            self.compact(v) if isinstance(v, Obj) else [self.compact(x.value) for x in v.items if isinstance(x.value, Obj)]
            if not getattr(v, 'new', False) or not v.items:
                continue
            if any(isinstance(x.value, Container) or getattr(x.value, 'inline', False) or '\n' in x.value.render()
                   or x.tail.strip() for x in v.items):
                continue
            inline = '{' + ', '.join(f'{x.keyraw}: {x.value.render()}' for x in v.items) + '}'
            if len(inline) > 72:
                continue
            ind = m.indent() or ''
            comments = [l.strip() for x in v.items for l in x.pre.split('\n') if l.strip()]
            if comments:
                m.pre = m.pre[:len(m.pre) - len(ind)] + ''.join(ind + c + '\n' for c in comments) + ind
            m.value = Scalar(inline)
            m.value.inline = True
            obj.dirty = True

    def switch(self, member, name):
        """0/1 -> false/true."""
        v = plain(member.value)
        if isinstance(v, bool):
            return member.value
        if v in (0, 1):
            self.note(f'{name}: {v} -> {"true" if v else "false"}')
            return Scalar('true' if v else 'false')
        raise ValueError(f'{name} must be true or false (or 0 or 1), not {v!r}')

    def method_shape(self, R):
        _, m = R.find('method')
        shape = None
        if m is not None:
            v = plain(m.value)
            low = str(v).lower()
            if low == 'cce':
                m.value = Scalar('"cCCE"')
                R.dirty = True
                self.note(f'method "{v}" -> "cCCE"')
            elif low in ('pcce', 'dsj', 'itb', 'dsjitb'):
                m.value = Scalar('"cCCE"')
                shape = low
                R.dirty = True
                self.note(f'method "{v}" -> "cCCE" with cluster_mode "{low}"')
        i, cm = R.find('cluster_mode')
        if cm is not None and str(plain(cm.value)).lower() == 'adaptive':
            cm.value = Scalar('"bundle"')
            R.dirty = True
            self.note('cluster_mode "adaptive" -> "bundle"')
        if shape is not None:
            if cm is not None:
                cm.value = Scalar(json.dumps(shape))
            else:
                i, _ = R.find('method')
                R.insert(i + 1, self.new_member('cluster_mode', json.dumps(shape), R.items[i].indent() or ''))
        self.shape = shape or (str(plain(cm.value)).lower() if cm is not None else 'standard')
        self.gcce = m is not None and str(plain(m.value)).lower() == 'gcce'

    def radii(self, R):
        for outer, inner in (('rbath', 'rbathcut'), ('rdip', 'rdipcut')):
            cut = R.get(inner)
            if cut is None:
                continue
            o = R.get(outer)
            v = plain(cut.value)
            if v == 0:
                R.drop(inner)
                self.note(f'{inner}: 0 removed (the default)')
                continue
            R.take(inner)
            if o is None:
                raise ValueError(f'{inner} without {outer}')
            o.value = Scalar(f'[{cut.value.render()}, {o.value.render()}]')
            if cut.tail.strip() and not o.tail.strip():
                o.tail = cut.tail
            self.note(f'{inner} -> {outer}: [inner, outer]')

    def renames(self, R):
        for old, new in (('bathadjust', 'bath_shift'), ('statefile', 'bath_state_file'),
                         ('nk', 'max_clusters_per_order'), ('addsubclus', 'add_subclusters')):
            _, m = R.find(old)
            if m is None:
                continue
            self.rename(m, new)
            if new == 'add_subclusters':
                m.value = self.switch(m, old)
            R.dirty = True
            self.note(f'{old} -> {new}')
        _, a = R.find('Adaptive')
        if a is not None and a.key != 'adaptive':
            self.rename(a, 'adaptive')
            R.dirty = True
            self.note('Adaptive -> adaptive')
        for name in ('DefectTotSpin', 'CorrTotSpin'):
            m = R.drop(name)
            if m is None:
                continue
            v = plain(m.value)
            if name == 'CorrTotSpin' and v != 0:
                raise ValueError('CorrTotSpin is not 0: it served old-format VASP tensor files, which have to be remade '
                                 'with example/tools/SpinTensorReader (format v3); the input cannot be converted as it is')
            self.note(f'{name} removed ({v}; ' + ('the qubit spin is used)' if name == 'DefectTotSpin' else 'it was 0)'))

    SAVEMODES = {'normal': 'coherence', 'info': 'coherence_tagged', 'all': 'cluster_coherence',
                 'allfull': 'cluster_density_matrix'}

    def savemode(self, R, top):
        """The savemode names before 2026-10-13 (normal, info, all, allfull) and the top-level save_nodiv, which CCEX
        refuses with their new form: "savemode": "<mode>", or {"mode": ..., "nodiv": ...} when save_nodiv was given."""
        i, m = R.find('savemode')
        if m is not None and not isinstance(m.value, Obj):
            v = plain(m.value)
            if isinstance(v, str) and v.lower() in self.SAVEMODES:
                m.value = Scalar('"' + self.SAVEMODES[v.lower()] + '"')
                R.dirty = True
                self.note(f'savemode "{v}" -> "{self.SAVEMODES[v.lower()]}"')
        j, n = R.find('save_nodiv')
        if n is None:
            return
        if m is not None and isinstance(m.value, Obj):
            raise ValueError('save_nodiv beside a savemode block: write nodiv inside the block')
        if not isinstance(plain(n.value), bool):     # 0/1 too: CCEX took only a JSON boolean, before and after A47
            raise ValueError('save_nodiv must be true or false')
        mode = m.value.render() if m is not None else '"coherence"'
        pos = i if i is not None else j
        if m is not None:
            R.take('savemode')
        R.take('save_nodiv')
        block = self.new_block(top)
        self.put(block, top, 'mode', value=mode)
        self.put(block, top, 'nodiv', value=n.value.render())
        mem = self.new_member('savemode', block, top, moved=m) if m is not None else self.new_member('savemode', block, top)
        mem.value, mem.shift = block, 0
        R.insert(min(pos, len(R.items)), mem)
        R.dirty = True
        self.note('save_nodiv -> savemode {mode, nodiv}')

    def move(self, R, top, block_key, pairs, position_keys=None, values=None):
        """Move R[old] to R[block_key][new] for (old, new) in pairs; returns the block (None if nothing moved)."""
        found = [(R.find(old)[0], old, new) for old, new in pairs if R.find(old)[1] is not None]
        if not found:
            return None
        pos = min(i for i, _, _ in found)
        first_pre = R.items[pos].pre                        # a block in place of a key written inline stays inline
        moved = []
        for _, old, new in sorted(found, key=lambda x: [p[0] for p in pairs].index(x[1])):
            _, m = R.take(old)
            moved.append((old, new, m))
        pos = min(pos, len(R.items))
        made = R.get(block_key) is None
        block = self.child_block(R, top, block_key, pos)
        if made and '\n' not in first_pre:
            R.items[pos].pre = first_pre
        for old, new, m in moved:
            value = values(old, m) if values else None
            if isinstance(new, tuple):                      # a block inside the block
                sub = self.child_block(block, self.indent_of(block, top), new[0])
                self.put(sub, self.indent_of(block, top), new[1], m, value, top)
                self.note(f'{old} -> {block_key}.{new[0]}.{new[1]}')
            else:
                self.put(block, top, new, m, value, top)
                self.note(f'{old} -> {block_key}.{new}')
        return block

    def blocks(self, R, top):
        pcce = [('sK', 'partition_size'), ('max_trial', 'trials'), ('max_iter', 'max_iterations'),
                ('kmeans_pp', 'kmeans_pp'), ('iter_detail', 'print_trials'), ('pcce_detail_limit', 'print_limit'),
                ('pcce_partition', 'partition'), ('pcce_connections', 'connections'),
                ('pcce_cache_mode', ('cache', 'mode')), ('pcce_cache_file', ('cache', 'file'))]
        if self.shape != 'pcce':
            for old, _ in pcce:
                m = R.get(old)
                if m is None or (old == 'pcce_cache_mode' and str(plain(m.value)).lower() != 'off'):
                    continue
                R.drop(old)
                self.note(f'{old} removed (it is read with cluster_mode pcce only)')
        self.move(R, top, 'pcce', pcce,
                  values=lambda old, m: self.switch(m, old) if old in ('kmeans_pp', 'iter_detail') else None)
        self.move(R, top, 'bundle', [('bundle_max_dim', 'max_dim'), ('bundle_merge_cutoff', 'merge_cutoff'),
                                     ('bundle_connections', 'connections'),
                                     ('bundle_connection_budget', 'connection_budget'),
                                     ('expansion_max_clusters', 'max_clusters'), ('expansion_max_work', 'max_work')])
        if R.get('adaptive_work_file') is not None:
            self.move(R, top, 'adaptive', [('adaptive_work_file', 'work_file')])
        self.move(R, top, 'pulse', [('npulse', 'count'), ('pulsename', 'name'), ('sequence', 'sequence'),
                                    ('pulseiter', 'per_qubit'), ('pulse_time', 'duration'),
                                    ('detuning_factor', 'detuning_factor')],
                  values=lambda old, m: self.switch(m, old) if old == 'pulseiter' else None)
        for name, inner, sub in (('gcce_refinement', 'gcce_refinement_max_sweeps', 'max_sweeps'),
                                 ('cluster_contributions', 'cluster_contribution_max_mb', 'max_mb')):
            i, m = R.find(name)
            _, n = R.find(inner)
            if m is not None and isinstance(m.value, Obj):
                continue
            if m is None and n is None:
                continue
            on = m.value.render() if m is not None else 'false'
            if m is not None and not isinstance(plain(m.value), bool):   # 0/1 too, as CCEX (readStrictBool)
                raise ValueError(f'{name} must be true or false')
            pos = i if i is not None else R.find(inner)[0]
            if m is not None:
                R.take(name)
            block = self.new_block(top)
            self.put(block, top, 'enabled', value=on)
            if n is not None:
                R.take(inner)
                self.put(block, top, sub, n, None, top)
            mem = self.new_member(name, block, top, moved=m) if m is not None else self.new_member(name, block, top)
            mem.value, mem.shift = block, 0
            R.insert(min(pos, len(R.items)), mem)
            R.dirty = True
            self.note(f'{name}' + (f' and {inner}' if n is not None else '') + f' -> {name} {{enabled, {sub}}}')

    def qubit(self, R, top):
        i, qf = R.find('qubitfile')
        if qf is not None:
            if R.get('Qubit') is not None:
                raise ValueError('the input has both qubitfile and a Qubit block (CCEX refused it before too)')
            self.move(R, top, 'Qubit', [('qubitfile', 'file'), ('qspin', 'spin'), ('qzfs', 'zfs'),
                                        ('qalphams', 'alphams'), ('qbetams', 'betams')])
            return
        _, q = R.find('Qubit')
        if q is None or not isinstance(q.value, Obj):
            for k in ('qspin', 'qzfs', 'qalphams', 'qbetams'):
                if R.get(k) is not None:
                    R.drop(k)
                    self.note(f'{k} removed (it was read with qubitfile only)')
            return
        Q = q.value
        qind = self.indent_of(Q, top)
        n = Q.drop('nqubit')
        if n is not None:
            qubits = Q.get('qubit')
            count = len(plain(qubits.value)) if qubits is not None else None
            if count != plain(n.value):
                raise ValueError(f'Qubit.nqubit is {plain(n.value)} but Qubit.qubit has {count} entries')
            self.note('Qubit.nqubit removed (the length of Qubit.qubit)')
        # the states: psia/psib were used when both were given, alphaidx/betaidx when both were given, else neither
        names = {'psia': None, 'psib': None, 'alphaidx': None, 'betaidx': None}
        first = None
        for k in names:
            j, m = Q.find(k)
            if m is not None:
                first = j if first is None else min(first, j)
        for k in names:
            names[k] = Q.take(k)[1]
        if first is not None:
            first = min(first, len(Q.items))
            if names['psia'] is not None and names['psib'] is not None:
                use = (('psia', 'alpha_state', 'vector'), ('psib', 'beta_state', 'vector'))
                dropped = [k for k in ('alphaidx', 'betaidx') if names[k] is not None]
            elif names['alphaidx'] is not None and names['betaidx'] is not None:
                use = (('alphaidx', 'alpha_state', 'index'), ('betaidx', 'beta_state', 'index'))
                dropped = [k for k in ('psia', 'psib') if names[k] is not None]
            else:
                use = ()
                dropped = [k for k in names if names[k] is not None]
            for k in dropped:
                self.note(f'Qubit.{k} removed (it was not used: the run took ' +
                          ('psia/psib' if use and use[0][2] == 'vector' else
                           'alphaidx/betaidx' if use else 'alphams/betams') + ')')
            for old, new, form in reversed(use):
                m = names[old]
                inner = self.new_block(qind)
                self.put(inner, qind, form, m, None, qind)
                mem = self.new_member(new, inner, qind)
                Q.insert(first, mem)
                Q.dirty = True
            for old, new, form in use:
                self.note(f'Qubit.{old} -> Qubit.{new}.{form}')
        self.qubit_overhaus = Q.take('overhaus')[1]

    def hamiltonian(self, R, top):
        common, ccce = [], []
        positions = []

        def take(key):
            i, m = R.find(key)
            if m is not None:
                positions.append(i)
                R.take(key)
            return m

        ovh = getattr(self, 'qubit_overhaus', None)
        if ovh is not None:
            common.append(('qubit_bath_mean_field', ovh, self.switch(ovh, 'Qubit.overhaus'), 'Qubit.overhaus'))
        vmf = take('vector_mean_field')
        if vmf is not None:
            v = plain(vmf.value)
            common.append(('bath_bath_mean_field', vmf, '"xyz"' if v in (True, 1) else '"z"', 'vector_mean_field'))
        m = R.get('hf_readmode')
        mode = plain(m.value) if m is not None else 0
        if mode not in (0, 1, 2, 3):
            raise ValueError(f'hf_readmode {mode} is not 0, 1, 2 or 3')
        if mode == 0:
            for k in ('hf_tensorfile', 'hf_ignore_oor'):
                if R.drop(k) is not None:
                    self.note(f'{k} removed (it is not read with hf_readmode 0)')
            if R.drop('hf_readmode') is not None:
                self.note('hf_readmode 0 removed (the default: hyperfine dipole "point", contact false)')
        else:
            hf = {k: take(k) for k in ('hf_readmode', 'hf_tensorfile', 'hf_ignore_oor')}
            block = self.new_block('')
            common.append(('hyperfine', hf['hf_readmode'], block, 'hf_readmode'))
            self.hf_block = (block, mode, hf)
        m = R.get('qd_readmode')
        mode = plain(m.value) if m is not None else 0
        if mode not in (0, 2, 3, 4):
            raise ValueError(f'qd_readmode {mode} is not 0, 2, 3 or 4')
        unread = ('qd_tensorfile', 'qd_tensorfile_woqubit', 'qd_cellpara') if mode == 0 else \
                 ('qd_tensorfile_woqubit', 'qd_cellpara') if mode == 2 else ()
        for k in unread:
            if R.drop(k) is not None:
                self.note(f'{k} removed (it is not read with qd_readmode {mode})')
        if mode == 0:
            if R.drop('qd_readmode') is not None:
                self.note('qd_readmode 0 removed (the default: quadrupole "off")')
        else:
            qd = {k: take(k) for k in ('qd_readmode', 'qd_tensorfile', 'qd_tensorfile_woqubit', 'qd_cellpara')}
            block = self.new_block('')
            common.append(('quadrupole', qd['qd_readmode'], block, 'qd_readmode'))
            self.qd_block = (block, mode, qd)
        for k in ('knight', 'hfmedi'):
            if self.gcce:
                if R.drop(k) is not None:
                    self.note(f'{k} removed (it has no effect with method gCCE)')
                continue
            m = take(k)
            if m is not None:
                ccce.append((k, m, self.switch(m, k), k))
        if not common and not ccce:
            return
        pos = min(positions) if positions else len(R.items)
        H = self.child_block(R, top, 'hamiltonian', min(pos, len(R.items)))
        hind = self.indent_of(H, top)
        for name, items in (('common', common), ('cCCE', ccce)):
            if not items:
                continue
            B = self.child_block(H, hind, name)
            bind = self.indent_of(B, hind)
            for key, moved, value, old in items:
                if isinstance(value, Obj):
                    value.close = '\n' + bind
                    self.fill(value, bind, key, top)
                self.put(B, hind, key, moved, value, top)
                self.note(f'{old} -> hamiltonian.{name}.{key}')

    def fill(self, block, bind, key, top):
        if key == 'hyperfine':
            _, mode, hf = self.hf_block
            self.put(block, bind, 'dipole', value='"dft"' if mode in (2, 3) else '"point"')
            self.put(block, bind, 'contact', value='true' if mode in (1, 3) else 'false')
            if mode != 0:
                if hf['hf_tensorfile'] is not None:
                    self.put(block, bind, 'file', hf['hf_tensorfile'], None, top)
                    self.note('hf_tensorfile -> hamiltonian.common.hyperfine.file')
                if hf['hf_ignore_oor'] is not None:
                    v = plain(hf['hf_ignore_oor'].value)
                    if v not in (0, 1):
                        raise ValueError(f'hf_ignore_oor {v} is not 0 or 1')
                    self.put(block, bind, 'missing_row', hf['hf_ignore_oor'], '"point_dipole"' if v else '"stop"', top)
                    self.note('hf_ignore_oor -> hamiltonian.common.hyperfine.missing_row')
        else:
            _, mode, qd = self.qd_block
            self.put(block, bind, 'structure', value=json.dumps({2: 'plain', 3: 'strained', 4: 'bubble'}[mode]))
            for old, new in (('qd_tensorfile', 'file'), ('qd_tensorfile_woqubit', 'file_without_defect'),
                             ('qd_cellpara', 'cell')):
                if qd[old] is not None and (mode != 2 or old == 'qd_tensorfile'):
                    self.put(block, bind, new, qd[old], None, top)
                    self.note(f'{old} -> hamiltonian.common.quadrupole.{new}')

    def defect(self, R, top):
        i, d = R.find('Defect')
        ex, ax = R.get('exstatefile'), R.get('avaaxfile')
        if d is None:
            for k in ('exstatefile', 'avaaxfile'):
                if R.get(k) is not None:
                    R.drop(k)
                    self.note(f'{k} removed (it is read with Defect species only)')
            return
        species = d.value
        if isinstance(species, Arr):
            renames = {'dfname': 'name', 'naddspin': 'n_nuclear_spins', 'navaax': 'n_configurations',
                       'eqs': 'quadrupole_moments', 'rxyzs': 'nuclear_positions', 'hypf': 'hyperfine'}
            for n, item in enumerate(species.items):
                e = item.value
                if not isinstance(e, Obj):
                    continue
                for m in list(e.items):
                    if m.key in renames:
                        self.rename(m, renames[m.key])
                        e.dirty = True
                if e.drop('apprx') is not None:
                    self.note(f'Defect[{n}].apprx removed (it changed nothing)')
            if any(isinstance(it.value, Obj) and it.value.dirty for it in species.items):
                self.note('Defect entries: dfname -> name, naddspin -> n_nuclear_spins, navaax -> n_configurations, '
                          'eqs -> quadrupole_moments, rxyzs -> nuclear_positions, hypf -> hyperfine')
            block = self.new_block(top)
            for key, old in (('state_file', 'exstatefile'), ('axis_file', 'avaaxfile')):
                if R.get(old) is not None:
                    _, m = R.take(old)
                    self.put(block, top, key, m, None, top)
                    self.note(f'{old} -> Defect.{key}')
            ind = self.indent_of(block, top)
            mem = Member('\n' + ind, '"species"', self.aligned(self.mid, self.width, '"species"'), species, '',
                         len(ind) - len(top))
            block.items.append(mem)
            d.value = block
            R.dirty = True
            self.note('Defect [...] -> Defect {"species": [...]}')
        elif isinstance(species, Obj):
            for key, old in (('state_file', 'exstatefile'), ('axis_file', 'avaaxfile')):
                if R.get(old) is not None:
                    _, m = R.take(old)
                    self.put(species, top, key, m, None, top)
                    self.note(f'{old} -> Defect.{key}')


def convert_text(text):
    c = Converter(text)
    new = c.convert()
    json.loads(strip_comments(new))          # what CCEX will read has to be JSON
    return new, c.log


def main(argv):
    mode = 'write'
    if argv and argv[0] in ('--check', '--stdout'):
        mode, argv = argv[0][2:], argv[1:]
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        return 0 if argv else 2
    status = 0
    for path in argv:
        try:                   # a file that cannot be read (missing, not UTF-8, ...) is reported like one that cannot be converted
            with open(path, newline='') as f:
                raw = f.read()
            crlf = '\r\n' in raw          # converted with \n and written back with \r\n (lines added by the converter too)
            text = raw.replace('\r\n', '\n') if crlf else raw
            new, log = convert_text(text)
        except Exception as error:
            print(f'{path}: not converted: {error}', file=sys.stderr)
            status = 2
            continue
        out = new.replace('\n', '\r\n') if crlf else new
        if mode == 'stdout':
            sys.stdout.write(out)
            continue
        if new == text:
            print(f'{path}: already in the new form')
            continue
        print(f'{path}:')
        for line in log:
            print(f'    {line}')
        if mode == 'check':
            status = max(status, 1)
            continue
        try:
            write_atomically(path, out)
        except Exception as error:
            print(f'{path}: not written (the file is unchanged): {error}', file=sys.stderr)
            status = 2
    return status


def write_atomically(path, text):
    """A temporary file beside the target, then os.replace: an interrupted write leaves the original whole."""
    target = os.path.realpath(path)       # through a symbolic link, as open(path, 'w') did
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(target), prefix='.' + os.path.basename(target) + '.')
    try:
        with os.fdopen(fd, 'w', newline='') as f:
            f.write(text)
        shutil.copymode(target, tmp)
        os.replace(tmp, target)
    except BaseException:
        os.unlink(tmp)
        raise


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
