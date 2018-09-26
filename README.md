# PhD Dissertation: The Bicycle Wheel
__Author:__ Matthew Ford

__Institution:__ Northwestern University

__Department:__ Mechanical Engineering

_Copyright 2018 by Matthew Ford - All Rights Reserved_

## How to build the thesis document (thesis.pdf)
1. Install a LaTeX distribution such as [TeXLive](https://www.tug.org/texlive/)
2. Install [Inkscape](https://inkscape.org) and make sure that the `inkscape` command is available from the command line
3. Download or clone this repository to your local computer
4. Open your favorite command line in the same directory as `thesis.tex`
5. Run `latexmk -pdf -quiet --shell-escape thesis.tex`

Alternatively, use your favorite LaTeX distribution and build with bibtex support. Make sure that your build system allows the latex command `\write18`, which is required to convert the figures from .svg to .pdf_tex.
