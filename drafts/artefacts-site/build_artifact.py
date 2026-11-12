"""Package the held-back read-artefact page as a claude.ai Artifact (internal sharing only).

Writes OUT/index.html (page body, no doctype: the Artifact host adds the skeleton), OUT/page.css
(site + page styles, dark mode also driven by the viewer's data-theme), OUT/artefacts.js and OUT/data/.
usage: build_artifact.py OUT
"""
import os, re, sys, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = sys.argv[1]
SITE = "https://kirk3gaard.github.io/MicroBench/"
os.makedirs(os.path.join(OUT, "data"), exist_ok=True)


def blocks(css):
    """Yield (prelude, body) for top-level rules, handling one level of nesting."""
    i, n = 0, len(css)
    while i < n:
        j = css.find("{", i)
        if j < 0: break
        depth, k = 1, j + 1
        while depth:
            depth += {"{": 1, "}": -1}.get(css[k], 0); k += 1
        yield css[i:j].strip(), css[j + 1:k - 1]
        i = k


def scope_dark(body):
    out = []
    for guard in (':root:not([data-theme="light"])', ':root[data-theme="dark"]'):
        for sel, decl in blocks(body):
            sels = []
            for s in sel.split(","):
                s = s.strip()
                sels.append(guard if s == ":root" else f"{guard} {s}")
            extra = " color-scheme: dark;" if sel.strip() == ":root" else ""
            rule = f'{", ".join(sels)} {{{decl}{extra}}}'
            out.append(rule if guard.startswith(":root[") else rule)
    media = [r for r in out[: len(out) // 2]]
    explicit = out[len(out) // 2:]
    return "@media (prefers-color-scheme: dark) {\n" + "\n".join(media) + "\n}\n" + "\n".join(explicit) + "\n"


def convert(css):
    css = re.sub(r"/\*.*?\*/", "", css, flags=re.S)
    res = []
    for pre, body in blocks(css):
        if pre.startswith("@media") and "prefers-color-scheme: dark" in pre:
            res.append(scope_dark(body))
        else:
            res.append(f"{pre} {{{body}}}")
    return "\n".join(res)


css = convert(open(os.path.join(ROOT, "docs", "style.css")).read()) + "\n" + convert(open(os.path.join(HERE, "artefacts.css")).read())
css += """
.draftbar { background: var(--draft-bg); color: var(--draft-fg); border-bottom: 1px solid var(--border); }
.draftbar .wrap { padding-block: 10px; font-size: .9rem; }
:root { --draft-bg: #fff4d6; --draft-fg: #5c4100; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --draft-bg: #3a2e0a; --draft-fg: #f5d77a; } }
:root[data-theme="dark"] { --draft-bg: #3a2e0a; --draft-fg: #f5d77a; }
"""
open(os.path.join(OUT, "page.css"), "w").write(css)

html = open(os.path.join(HERE, "artefacts.html")).read()
body = html.split("<body>", 1)[1].rsplit("</body>", 1)[0]
body = body.replace('href="./"', f'href="{SITE}"').replace('href="analyses/"', f'href="{SITE}analyses/"')
draft = ('<div class="draftbar"><div class="wrap"><b>Internal draft, not yet verified.</b> These numbers come from an '
         'automated scan of every read and have not been checked by anyone else yet. Please do not share outside the group.</div></div>\n')
head = ('<title>MicroBench Read Artefacts</title>\n'
        '<link rel="icon" href="microbench.jpg">\n<link rel="stylesheet" href="page.css">\n')
open(os.path.join(OUT, "index.html"), "w").write(head + draft + body.strip() + "\n")
shutil.copy(os.path.join(HERE, "artefacts.js"), OUT)
shutil.copy(os.path.join(ROOT, "docs", "microbench.jpg"), OUT)
for f in os.listdir(os.path.join(HERE, "data")):
    shutil.copy(os.path.join(HERE, "data", f), os.path.join(OUT, "data", f))
print(sorted(os.listdir(OUT)), sorted(os.listdir(os.path.join(OUT, "data"))))
