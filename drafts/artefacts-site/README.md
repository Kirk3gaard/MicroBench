# Read-artefact page (held back)

The interactive read-artefact page is kept out of `docs/` (and so off GitHub Pages) until its
findings have been checked further. It is the same page that was briefly live at
`artefacts.html`; data for ZymoHMW, ZymoGut and ZymoOral are in `data/`.

Preview locally (the page uses the site's `style.css` and logo):

    mkdir -p /tmp/preview && cp -r docs/* drafts/artefacts-site/* /tmp/preview/
    cd /tmp/preview && python3 -m http.server 8000   # open http://localhost:8000/artefacts.html

To publish it, move the files back: `artefacts.html/.js/.css` to `docs/`, `data/*` to `docs/data/`,
restore the link in `docs/index.html`, and the `EXTRA` entry and footer link in
`scripts/build_analysis_pages.py`. New runs: `analysis/artefacts/merge_page.py ... drafts/artefacts-site/data`.
