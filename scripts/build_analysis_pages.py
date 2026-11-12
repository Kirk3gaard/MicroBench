#!/usr/bin/env python3
"""Build web pages for the analyses in analysis/*/ from their knitted README.md.

For every analysis folder with a README.md (the github_document output of README.Rmd) this writes
docs/analyses/<name>.html, copies the figures to docs/analyses/<name>/, and writes an overview at
docs/analyses/index.html. Rerun after re-knitting a README; the GitHub Action does it on push.

Requires: pip install markdown
"""
import os, re, shutil, html, datetime
import markdown

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "analysis")
OUT = os.path.join(ROOT, "docs", "analyses")
REPO = "https://github.com/Kirk3gaard/MicroBench"

# extra entries that live elsewhere on the site
EXTRA = []  # the read-artefact page is held back in drafts/artefacts-site until verified
# nicer titles/summaries for knitted reports (fallback: first heading and the Aim section)
META = {"zymohmw": {"title": "Basecalling models and assemblers on the ZymoHMW mock",
                    "tags": ["basecalling", "accuracy", "assembly"]}}

SLUG = lambda s: re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")

def parse_readme(path):
    text = open(path, encoding="utf-8").read()
    # github_document header: "README\n====\nAuthor\nDate"
    m = re.match(r"^(.+)\n=+\n(.*?)\n\n", text, re.S)
    author, date = "", ""
    if m:
        lines = [x.strip() for x in m.group(2).splitlines() if x.strip()]
        author = lines[0] if lines else ""
        date = lines[1] if len(lines) > 1 else ""
        text = text[m.end():]
    text = re.sub(r"<!-- -->", "", text)
    title = re.search(r"^# (.+)$", text, re.M)
    title = title.group(1).strip() if title else os.path.basename(os.path.dirname(path))
    text = re.sub(r"^# .+\n", "", text, count=1, flags=re.M)
    aim = re.search(r"^## Aim\s*\n+(.+?)(?:\n## |\Z)", text, re.S | re.M)
    concl = re.search(r"^## Conclusion\s*\n+(.+?)(?:\n## |\Z)", text, re.S | re.M)
    return {"title": title, "author": author, "date": date, "body": text,
            "aim": " ".join(aim.group(1).split()) if aim else "",
            "conclusion": " ".join(concl.group(1).split()) if concl else ""}

def render_body(md, name):
    figs = []
    def img(m):
        src = m.group(1)
        figs.append(src)
        return f"![]({name}/{os.path.basename(src)})"
    md = re.sub(r"!\[[^\]]*\]\(([^)]+\.(?:png|jpg|jpeg|svg|gif))\)", img, md)
    body = markdown.markdown(md, extensions=["tables", "fenced_code", "sane_lists"])
    # headings get ids; figures become cards titled by the nearest heading above
    toc, current = [], {"h2": "", "h3": ""}
    used = set()
    def heading(m):
        lvl, inner = m.group(1), m.group(2)
        plain = re.sub(r"<[^>]+>", "", inner)
        hid = SLUG(plain) or "section"
        base, k = hid, 2
        while hid in used: hid = f"{base}-{k}"; k += 1
        used.add(hid)
        if lvl in ("2", "3"): toc.append((int(lvl), plain, hid))
        return f'<h{lvl} id="{hid}">{inner}<a class="anchor" href="#{hid}" aria-label="Link to this section">#</a></h{lvl}>'
    body = re.sub(r"<h([1-6])>(.*?)</h\1>", heading, body)
    caps = []
    for part in re.split(r'(<h[23] id="[^"]+">.*?</h[23]>)', body):
        clean = lambda p: re.sub(r"<[^>]+>", "", re.sub(r'<a class="anchor".*?</a>', "", p)).strip()
        if part.startswith("<h2"): current["h2"] = clean(part); current["h3"] = ""
        elif part.startswith("<h3"): current["h3"] = clean(part)
        else: caps.extend([" – ".join(x for x in (current["h2"], current["h3"]) if x)] * part.count("<img "))
    it = iter(caps)
    def figure(m):
        cap = html.escape(next(it, ""))
        src = m.group(1)
        return (f'<figure class="fig"><a href="{src}" class="zoom" data-cap="{cap}"><img src="{src}" alt="{cap}" loading="lazy"></a>'
                f'<figcaption>{cap} <a href="{src}" download>PNG</a></figcaption></figure>')
    body = re.sub(r'<p><img alt="" src="([^"]+)" ?/?></p>', figure, body)
    # external links open in a new tab
    body = re.sub(r'<a href="(https?://[^"]+)"', r'<a href="\1" target="_blank" rel="noopener"', body)
    return body, toc, figs

PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title} – MicroBench</title>
<meta name="description" content="{desc}">
<link rel="icon" href="{up}microbench.jpg">
<link rel="stylesheet" href="{up}style.css">
<link rel="stylesheet" href="{up}analyses.css">
</head>
<body>
<header class="hero slim">
  <div class="wrap hero-inner">
    <a href="{up}./"><img src="{up}microbench.jpg" alt="MicroBench home" class="logo" width="64" height="64"></a>
    <div>
      <nav class="crumbs"><a href="{up}./">Datasets</a> <span>›</span> {crumb}</nav>
      <h1>{title}</h1>
      {sub}
    </div>
  </div>
</header>
<main class="wrap">
{main}
</main>
<footer class="wrap muted">MicroBench · Aalborg University · <a href="{up}./">Datasets</a> · <a href="{up}analyses/">Analyses</a></footer>
<div class="lightbox" id="lb" hidden><button type="button" aria-label="Close">×</button><img alt=""><p></p></div>
<script>
const lb = document.getElementById("lb");
document.addEventListener("click", (e) => {{
  const a = e.target.closest("a.zoom");
  if (a) {{ e.preventDefault(); lb.querySelector("img").src = a.href; lb.querySelector("p").textContent = a.dataset.cap; lb.hidden = false; return; }}
  if (!lb.hidden && (e.target === lb || e.target.closest("#lb button"))) lb.hidden = true;
}});
document.addEventListener("keydown", (e) => {{ if (e.key === "Escape") lb.hidden = true; }});
const tw = document.getElementById("tocwrap");
if (tw && window.matchMedia("(min-width: 861px)").matches) tw.open = true;
const links = [...document.querySelectorAll(".toc a")];
if (links.length && "IntersectionObserver" in window) {{
  const io = new IntersectionObserver((es) => es.forEach((en) => {{
    if (en.isIntersecting) links.forEach((l) => l.classList.toggle("on", l.getAttribute("href") === "#" + en.target.id));
  }}), {{ rootMargin: "0px 0px -70% 0px" }});
  document.querySelectorAll("article h2[id], article h3[id]").forEach((h) => io.observe(h));
}}
</script>
</body>
</html>
"""

def build():
    os.makedirs(OUT, exist_ok=True)
    entries = []
    for name in sorted(os.listdir(SRC)):
        d = os.path.join(SRC, name)
        rd = os.path.join(d, "README.md")
        if not os.path.isfile(rd) or name == "artefacts":
            continue
        info = parse_readme(rd)
        meta = META.get(name, {})
        title = meta.get("title", info["title"])
        body, toc, figs = render_body(info["body"], name)
        fdir = os.path.join(OUT, name)
        if os.path.isdir(fdir): shutil.rmtree(fdir)
        os.makedirs(fdir)
        for f in figs:
            shutil.copy(os.path.join(d, f), os.path.join(fdir, os.path.basename(f)))
        files = sorted(x for x in os.listdir(d) if os.path.isfile(os.path.join(d, x)) and not x.startswith("."))
        code_links = "".join(f'<li><a href="{REPO}/blob/main/analysis/{name}/{html.escape(x)}" target="_blank" rel="noopener">{html.escape(x)}</a></li>'
                             for x in files if x in ("Snakefile", "config.yaml", "README.Rmd", "environment.yml", "create_env.sh") or x.endswith(".sh"))
        tochtml = "".join(f'<li class="l{l}"><a href="#{i}">{html.escape(t)}</a></li>' for l, t, i in toc)
        key = ""
        main = f"""<div class="report">
  <aside class="side">
    <details class="tocwrap" id="tocwrap"><summary>Contents</summary><nav class="toc" aria-label="Contents"><ul>{tochtml}</ul></nav></details>
    <div class="sidebox"><b>Code</b><ul>{code_links}</ul>
      <p><a href="{REPO}/tree/main/analysis/{name}" target="_blank" rel="noopener">Folder on GitHub ↗</a></p></div>
  </aside>
  <article class="prose">{key}{body}</article>
</div>"""
        sub = f'<p class="lede">{html.escape(info["author"])}{" · updated " + html.escape(info["date"]) if info["date"] else ""} · {len(figs)} figures</p>'
        page = PAGE.format(title=html.escape(title), desc=html.escape(info["aim"] or title), up="../",
                           crumb=f'<a href="./">Analyses</a> <span>›</span> <b>{html.escape(name)}</b>', sub=sub, main=main)
        open(os.path.join(OUT, f"{name}.html"), "w", encoding="utf-8").write(page)
        entries.append({"slug": f"{name}.html", "title": title, "summary": info["aim"] + (" " + info["conclusion"] if info["conclusion"] else ""),
                        "code": f"analysis/{name}", "tags": meta.get("tags", []), "date": info["date"], "thumb": f"{name}/{os.path.basename(figs[0])}" if figs else None,
                        "figs": len(figs)})
    entries += EXTRA
    cards = []
    for e in entries:
        thumb = f'<img src="{e["thumb"]}" alt="" loading="lazy">' if e.get("thumb") else '<div class="thumb-ph">interactive</div>'
        tags = "".join(f'<span class="tag">{html.escape(t)}</span>' for t in e.get("tags", []))
        meta = " · ".join(x for x in (e.get("date", ""), f'{e["figs"]} figures' if e.get("figs") else "") if x)
        cards.append(f"""<a class="acard" href="{e['slug']}">
  <div class="thumb">{thumb}</div>
  <div class="abody"><h2>{html.escape(e['title'])}</h2><p>{html.escape(e['summary'])}</p>
  <div class="ameta">{tags}<span class="muted">{html.escape(meta)}</span></div></div></a>""")
    main = f"""<p class="lede-in">Analyses run on the MicroBench data. Each page is built from the analysis folder in the repository, where the full workflow (Snakefile, environments, R Markdown) lives.</p>
<div class="agrid">{''.join(cards)}</div>"""
    page = PAGE.format(title="Analyses", desc="Analyses of the MicroBench nanopore benchmark data.", up="../",
                       crumb="<b>Analyses</b>", sub="", main=main)
    open(os.path.join(OUT, "index.html"), "w", encoding="utf-8").write(page)
    print("built", [e["slug"] for e in entries])

if __name__ == "__main__":
    build()
