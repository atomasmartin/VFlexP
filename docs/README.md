# Documentation

Two formats, two readers.

| path | for | start at |
|---|---|---|
| [`index.html`](index.html) | people: the interactive HTML documentation site, with figures, Simulink diagrams, and rendered equations | open [`index.html`](index.html) or [`doc_Home.html`](doc_Home.html) in a browser |
| [`agent/`](agent/) | AI agents, scripts, and terminals: the same material as plain Markdown, plus file formats, `VFlexPSession` / `core/` API contracts, UI architecture, and modular Simulink library specs | [`agent/README.md`](agent/README.md) |

The HTML site needs no build step and no server: every page (`doc_*.html`) is a plain HTML file whose sidebar navigation, table of contents, and light/dark theme come from `scripts/load_sidebar.js`, `scripts/load_theme.js`, and `scripts/switch_themes.js`. Equations are rendered by MathJax from a CDN.

## Agent documentation (`docs/agent/`)

| file | what it holds |
|---|---|
| [`agent/README.md`](agent/README.md) | what VFlexP is, how to drive it interactively (`vflexp`) or headlessly (`VFlexPSession` + `core/`), workflow pipeline, `git diff` traceability for `.slx`, and repository map |
| [`agent/formats.md`](agent/formats.md) | the three-file case format: workbook (`.xlsx`), dynamic data (`.dyr`), case settings (`.json`), `$name` symbolic variables, and simulation disturbance events |
| [`agent/models.md`](agent/models.md) | equations, per-unit bases, reference frames, complexity vocabulary, and every modular Simulink device model in `Simulink/lib/` |
| [`agent/reference.md`](agent/reference.md) | complete reference of `VFlexPSession`, `core/vflexp_*.m`, `ui/`, `functions/`, `tools/`, and `Simulink/lib/` mask parameters |
