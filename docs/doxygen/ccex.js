// CCEX on top of Doxygen's navigation tree.
//
// 1. Doxygen expands the tree entry of the page being shown only when that page is
//    an index page; a page entry with children of its own ("Input file") stays
//    shut. And the root of the tree (the project name) is the main page too, so on
//    the main page Doxygen marks the root, not the "Overview" entry under it. Once
//    the tree has marked the current page: move the mark from the root to the
//    entry that shows the same page, then open the marked entry's children.
// 2. The sidebar folds to a narrow rail with a button, and remembers that across
//    pages; the rail shows the path of the page being read and unfolds on click.
// 3. The sidebar's dragged width goes into the theme's width variable, so that
//    the sidebar, the header, the search box and the page follow it together.
// 4. The table of contents beside a page opens, at every level, only the entry that
//    holds the section being read (the theme marks its link "active").
(function () {
  function fileOf(href) { return (href || '').split('#')[0].split('/').pop(); }
  function rootEntry() { return document.querySelector('#nav-tree-contents > ul > li'); }
  function labelOf(li) {
    // the entry's link sits in span.label inside div.item
    var link = li.querySelector(':scope > .item a[href]:not([href^="javascript"]), :scope > .item a.nolink');
    return link ? link.textContent.trim() : '';
  }

  // the rail's text: the marked entry and its ancestors, "Input file › Cluster tags"
  function showPath(item) {
    var rail = document.getElementById('ccex-nav-rail');
    if (!rail) return;
    var names = [], root = rootEntry();
    for (var li = item.parentNode; li && li !== root && li.tagName === 'LI'; li = li.parentNode.parentNode) {
      var name = labelOf(li);
      if (name) names.unshift(name);
      if (!li.parentNode || li.parentNode.tagName !== 'UL') break;
    }
    rail.firstChild.textContent = names.length ? names.join(' › ') : document.title.replace(/^CCEX:\s*/, '');
  }

  function openSelected() {
    var item = document.querySelector('#nav-tree .item.selected');
    if (!item) return false;
    var root = rootEntry();
    if (root && item.parentNode === root) {
      var here = fileOf(location.pathname) || 'index.html';
      var entries = root.querySelectorAll(':scope > ul > li > .item');
      for (var i = 0; i < entries.length; i++) {
        var link = entries[i].querySelector('a[href]:not([href^="javascript"])');
        if (link && fileOf(link.getAttribute('href')) === here) {
          item.classList.remove('selected');
          entries[i].classList.add('selected');
          item = entries[i];
          break;
        }
      }
    }
    var arrow = item.querySelector('span.arrow');
    if (arrow && arrow.textContent.trim() === '►') arrow.parentNode.click();
    showPath(item);
    return true;
  }

  var KEY = 'ccex-nav-hidden';
  function remembered() { try { return localStorage.getItem(KEY) === '1'; } catch (e) { return false; } }
  function remember(hidden) { try { localStorage.setItem(KEY, hidden ? '1' : '0'); } catch (e) {} }
  function showState(button) {
    var hidden = document.documentElement.classList.contains(KEY);
    button.textContent = hidden ? '▶' : '◀';
    button.title = hidden ? 'show the contents' : 'fold the contents';
  }
  if (remembered()) document.documentElement.classList.add(KEY);

  document.addEventListener('DOMContentLoaded', function () {
    var button = document.createElement('button');
    button.id = 'ccex-nav-toggle';
    button.type = 'button';
    button.addEventListener('click', function () {
      var hidden = document.documentElement.classList.toggle(KEY);
      remember(hidden);
      showState(button);
    });
    document.body.appendChild(button);
    showState(button);

    var rail = document.createElement('div');
    rail.id = 'ccex-nav-rail';
    rail.title = 'show the contents';
    rail.appendChild(document.createElement('span'));
    rail.firstChild.textContent = document.title.replace(/^CCEX:\s*/, '');
    rail.addEventListener('click', function () { button.click(); });
    document.body.appendChild(rail);

    // the theme keeps the sidebar at the width variable; a drag (jQuery UI's
    // resize event, from Doxygen's resize.js) moves the variable, and the width
    // is remembered for the next page
    var WIDTH = 'ccex-nav-width';
    function setWidth(w) { document.documentElement.style.setProperty('--side-nav-fixed-width', w + 'px'); }
    try { var stored = parseInt(localStorage.getItem(WIDTH), 10); if (stored >= 200) setWidth(stored); } catch (e) {}
    if (window.jQuery) {
      jQuery('#side-nav').on('resize', function (event, ui) {
        if (!ui || !ui.size) return;
        var w = Math.max(200, Math.min(Math.round(ui.size.width), Math.round(window.innerWidth * 0.7)));
        setWidth(w);
        try { localStorage.setItem(WIDTH, String(w)); } catch (e) {}
      });
    }

    var tries = 0;
    var timer = setInterval(function () {
      if (openSelected() || ++tries > 50) clearInterval(timer);
    }, 100);
  });
})();

// 4. (above) the table of contents opens the section being read and folds the others
(function () {
  // at every level, the entry that holds the section being read is open and its siblings are folded
  function sync(toc) {
    var lists = [toc.querySelector(':scope > ul')];
    while (lists.length) {
      var ul = lists.pop();
      if (!ul) continue;
      var items = ul.querySelectorAll(':scope > li');
      var open = null;
      for (var i = 0; i < items.length; i++) if (items[i].querySelector('a.active')) { open = items[i]; break; }
      if (!open && ul.parentNode === toc && items.length) open = items[0];
      for (var j = 0; j < items.length; j++) {
        items[j].classList.toggle('ccex-open', items[j] === open);
        var sub = items[j].querySelector(':scope > ul');
        if (sub) lists.push(sub);
      }
    }
  }
  document.addEventListener('DOMContentLoaded', function () {
    var toc = document.querySelector('.contents > .toc');
    if (!toc) return;
    sync(toc);
    new MutationObserver(function () { sync(toc); })
      .observe(toc, {subtree: true, attributes: true, attributeFilter: ['class']});
  });
})();
