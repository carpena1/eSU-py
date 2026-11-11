#!/bin/sh
# Rebuild eSU_Python_Guide.pdf from eSU_Python_Guide.md.
# Needs pandoc (https://pandoc.org) and Google Chrome (or Chromium).
# Usage (from the docs/ folder):  sh build_pdf.sh
set -e
cd "$(dirname "$0")"
CHROME="${CHROME:-/Applications/Google Chrome.app/Contents/MacOS/Google Chrome}"
pandoc eSU_Python_Guide.md -s --toc --toc-depth=2 --math-method=mathml \
       --css guide.css --embed-resources -o guide.html
PROFILE="$(mktemp -d)"                       # throwaway Chrome profile
"$CHROME" --headless=new --disable-gpu --no-pdf-header-footer \
          --user-data-dir="$PROFILE" --print-to-pdf="$PWD/eSU_Python_Guide.pdf" "file://$PWD/guide.html"
rm -rf guide.html "$PROFILE"
echo "wrote eSU_Python_Guide.pdf"
