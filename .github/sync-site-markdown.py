"""Keep public markdown companions in sync with generated HTML, without fits."""
import argparse
from html.parser import HTMLParser
from pathlib import Path
import subprocess


class MainContent(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.inside = False
        self.parts = []

    def handle_starttag(self, tag, attrs):
        if tag == "main":
            self.inside = True
        elif self.inside:
            self.parts.append(self.get_starttag_text())

    def handle_endtag(self, tag):
        if tag == "main":
            self.inside = False
        elif self.inside:
            self.parts.append(f"</{tag}>")

    def handle_data(self, data):
        if self.inside:
            self.parts.append(data)

    def handle_entityref(self, name):
        if self.inside:
            self.parts.append(f"&{name};")

    def handle_charref(self, name):
        if self.inside:
            self.parts.append(f"&#{name};")


def synchronize(root):
    # Old nested article aliases must not continue serving superseded claims.
    aliases = root / "articles/articles"
    alias_count = 0
    if aliases.exists():
        for path in aliases.glob("*.html"):
            target = path.parent.parent / path.name
            if target.is_file():
                url = "../" + path.name
                path.write_text(
                    '<!doctype html><html><head><meta charset="utf-8">'
                    f'<link rel="canonical" href="{url}">'
                    f'<meta http-equiv="refresh" content="0; url={url}">'
                    '</head><body><main><h1>Moved</h1>'
                    f'<p><a href="{url}">Read the maintained article.</a></p>'
                    '</main></body></html>\n'
                )
                alias_count += 1
    count = 0
    for path in sorted(root.rglob("*.html")):
        parser = MainContent()
        parser.feed(path.read_text())
        if not parser.parts:
            continue
        result = subprocess.run(
            ["pandoc", "--from=html", "--to=gfm", "--wrap=none"],
            input="".join(parser.parts), text=True, check=True, capture_output=True)
        # Preserve Markdown hard breaks without trailing whitespace in Git diffs.
        text = "\n".join(line.rstrip() + ("<br>" if line.endswith("  ") else "")
                         for line in result.stdout.splitlines()) + "\n"
        path.with_suffix(".md").write_text(text)
        count += 1
    print(f"Synchronized {count} markdown companions; redirected {alias_count} legacy article aliases.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("site", type=Path)
    synchronize(parser.parse_args().site)
