# docs/

The documentation, also published at <https://cqml-skku.github.io/CCEX/>.

| Folder | Contents |
|---|---|
| `user/` | the manual. Start with `Overview.md` and `GettingStarted.md` (building, running, the command line, the environment variables). `InputFile.md` and its pages describe every input key with its type, unit and default; `OutputFiles.md` describes every output file. `Examples.md` indexes the tutorials, `Analysis.md` and `Tools.md` cover the programs in `example/`, and `Changes.md` summarizes each version (`CHANGELOG.md` has the full list). |
| `methods/` | one page per method: cCCE and gCCE, pCCE, bundles and the adaptive refinement, the master-equation CCE, DEER bath pulses, the coordinate frame rotation, the Gaussian phase noise, the pseudospin model and the GPU backend. `Methods.md` lists them. |
| `doxygen/` | the website. `doxygen Doxyfile` (or `do_doxy`) builds it into `html/` from the pages above, the tutorial READMEs and the source. `mdfilter.sh` adapts the GitHub links and formulas for Doxygen; the look is the doxygen-awesome theme with `ccex.css` on top. |
| `images/` | the figures the pages use; `images/tutorials/` links each tutorial's `expected.png`. |

The checks are described in [test/README.md](../test/README.md) and the source layout in
[src/README.md](../src/README.md).
