%%% This is nuthesis.cls
%%% for writing doctoral dissertations at the 
%%% Mathematics Department of Northwestern University; and
%%% Electrical Engineering and Computer Science, Northwestern University
%%% 
%%% Author: Miguel A. Lerma (mlerma@math.northwestern.edu)
%%% Date:5/4/2005
%%% Edited by: Debjit Sinha (debjitsinha@yahoo.com) 
%%%	for NU Electronic dissertation requirements.
%%% Date:4/15/2006
%%%
%%% ====================================================================
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\NeedsTeXFormat{LaTeX2e}% LaTeX 2.09 can't be used (nor non-LaTeX)
[1994/12/01]% LaTeX date must December 1994 or later
\ProvidesClass{nuthesis}[2005/05/04 v1.6a by Miguel A. Lerma]
\PassOptionsToClass{oneside}{amsbook}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{amsbook}}
\ProcessOptions\relax
\LoadClass{amsbook}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%\def\singlespace{\renewcommand{\baselinestretch}{1}}
%\def\singlehalfspace{\renewcommand{\baselinestretch}{1.4}}
%\def\doublespace{\renewcommand{\baselinestretch}{1.8}}

%\RequirePackage{doublespace}       % MAL - 5/8/2001
%\def\onehalfspace{\spacing{1.4}}   % MAL - 5/8/2001
%\def\doublespace{\spacing{1.8}}    % MAL - 5/8/2001

\RequirePackage{setspace}                        % MAL - 4/19/2003
\renewcommand{\onehalfspacing}{\setstretch{1.4}} % MAL - 4/19/2003
\renewcommand{\doublespacing}{\setstretch{1.8}}  % MAL - 4/19/2003

% Page numbering should be at normal size 
% (amsbook.cls yields scriptsize) - MAL 6/17/2003
\def\ps@plain{\ps@empty
  \def\@oddfoot{\normalfont\normalsize \hfil\thepage\hfil}%
  \let\@evenfoot\@oddfoot}

% Another place where amsbook.cls yields scriptsize
% (amsbook.cls yields scriptsize) - MAL 7/15/2003
\def\ps@headings{\ps@empty
  \def\@evenhead{%
    \setTrue{runhead}%
    \normalfont\normalsize        % changed from scriptsize - MAL 7/15/2003
    \rlap{\thepage}\hfil \leftmark{}{}\hfil}%
  \def\@oddhead{%
    \setTrue{runhead}%
    \normalfont\normalsize \hfil  % changed from scriptsize - MAL 7/15/2003
    \rightmark{}{}\hfil \llap{\thepage}}%
  \let\@mkboth\markboth
  \def\partmark{\@secmark\markboth\partrunhead\partname}%
  \def\chaptermark{%
    \@secmark\markboth\chapterrunhead{}}%
  \def\sectionmark{%
    \@secmark\markright\sectionrunhead\sectionname}%
}

\def\contentsname{Table of Contents}

\def\chapter{%
  \if@openright\cleardoublepage\else\clearpage\fi
  \thispagestyle{myheadings}\global\@topnum\z@
  \@afterindenttrue \secdef\@chapter\@schapter}

\setlength{\topmargin}{0.0in}
\setlength{\headheight}{0.2in}
\setlength{\headsep}{0.5in}

\setlength{\textheight}{8.2in}
\setlength{\footskip}{0.0in}

\setlength{\textwidth}{6.5in}
\setlength{\oddsidemargin}{0.0in}
\setlength{\evensidemargin}{0.0in}

\renewcommand\thesection{\thechapter.\@arabic\c@section}
\renewcommand\thefigure {\thechapter.\arabic{figure}}
\renewcommand\thetable {\thechapter.\arabic{table}}

%% SECTION heading spacing and font
%%
% arguments are: #1 - sectiontype name
% (for \@sect)   #2 - section level
%                #3 - indent
%                #4 - top separation
%                #5 - tempskipa, bottom separation
%                #6 - font control
% if tempskipa <= 0 the heading is a run-in heading

% \def\section{\@startsection{section}{1}%
  % \z@{.7\linespacing\@plus\linespacing}{.5\linespacing}%
  % {\normalfont\bfseries\centering}}
\def\section{\@startsection{section}{1}%
  \z@{.7\linespacing\@plus\linespacing}{.5\linespacing}%
  {\normalfont\bfseries}}
\def\subsection{\@startsection{subsection}{2}%
  \z@{-3.25ex\@plus -1ex \@minus -.2ex}{1.5ex \@plus .2ex}%
  {\normalfont\bfseries}}
\def\subsubsection{\@startsection{subsubsection}{3}%
  \z@\z@{-.5em}%
  {\normalfont\bfseries}}


\def\th@plain{%
  \thm@headfont{\bfseries}
  \let\thmhead\thmhead@plain \let\swappedhead\swappedhead@plain
  \thm@preskip.5\baselineskip\@plus.2\baselineskip
                                    \@minus.2\baselineskip
  \thm@postskip\thm@preskip
  \itshape
}
\def\th@definition{%
  \thm@headfont{\bfseries}%
  \let\thmhead\thmhead@plain \let\swappedhead\swappedhead@plain
  \thm@preskip.5\baselineskip\@plus.2\baselineskip
                                    \@minus.2\baselineskip
  \thm@postskip\thm@preskip
  \upshape
}
\def\th@remark{%
  \thm@headfont{\bfseries}%
  \let\thmhead\thmhead@plain \let\swappedhead\swappedhead@plain
  \thm@preskip.5\baselineskip\@plus.2\baselineskip
                                    \@minus.2\baselineskip
  \thm@postskip\thm@preskip
  \upshape
}


% This redefinition of the proof environment is obsolete - MAL 5/4/2005
%\renewenvironment{proof}[1][\proofname]{\par \normalfont
%  \topsep6\p@\@plus6\p@ \trivlist \itemindent\z@
%  \item[\hskip\labelsep\bfseries
%    #1\@addpunct{.}]\ignorespaces
%}{%
%  \qed\endtrivlist
%}

% Better just redefine the font for proofname 
% - suggested by Ezra Getzler 5/4/2005
\renewcommand{\proofname}{\bfseries{Proof}}


\numberwithin{equation}{chapter}


\def\@captionheadfont{\normalfont}  % \scshape in amsbook.cls


\def\@title{Title of dissertation required}

\def\degree#1{\def\d@gree{#1}}
\def\d@gree{DOCTOR OF PHILOSOPHY}   % spelling corrected - MAL 3/8/2002
\def\field#1{\def\fi@ld{#1}}
\def\fi@ld{Mathematics}
\def\campus#1{\def\c@mpus{#1}}
\def\c@mpus{EVANSTON, ILLINOIS}

\def\graduation@month{\ifcase\month\or
    June\or June\or June\or June\or June\or June\or December\or
    December\or December\or December\or December\or December\fi}
\def\graduation@year{\number\year}

\def\graduationmonth#1{\def\graduation@month{#1}}
\def\graduationyear#1{\def\graduation@year{#1}}

\def\@date{\graduation@month \graduation@year}



\def\frontmatter{%
  \cleardoublepage
%  \pagenumbering{roman}
%  \pagestyle{plain}
  \pagenumbering{arabic}
  \pagestyle{myheadings}
}
\def\mainmatter{%
  \cleardoublepage
%  \pagenumbering{arabic}
%  \pagestyle{myheadings}
}
\def\backmatter{\cleardoublepage}


\def\titlepage{%
  \cleardoublepage
  \thispagestyle{empty}\setcounter{page}{1}}
\def\endtitlepage{\newpage}


\def\copyrightpage{%
  \cleardoublepage
\begin{center}
  \ \\
  \vspace{\fill}
\copyright \ Copyright by \authors{} \graduation@year \\
All Rights Reserved
\vspace{\fill}
\end{center}}
\def\endcopyrightpage{\newpage}


\def\maketitle{\par
  \@topnum\z@ % this prevents figures from falling at the top of page 1
\begin{titlepage}
\begin{center}
NORTHWESTERN UNIVERSITY\par
\vfill
{\large \@title}\par
\vfill
A DISSERTATION\par
\vspace{0.5cm}
SUBMITTED TO THE GRADUATE SCHOOL\par
IN PARTIAL FULFILLMENT OF THE REQUIREMENTS\par
\vspace{0.5cm}
for the degree\par
\vspace{0.5cm}
\d@gree{}\par
\vfill
Field of \fi@ld\par
\vfill
By\par
\vspace{0.5cm}
\ifx\@empty\authors
Author's name required\par
  \else
  \authors\par
  \fi
\vfill
\c@mpus\par
\vspace{0.5cm}
\graduation@month{} \graduation@year{}
\end{center}
\end{titlepage}
}


\renewcommand{\abstract}{%
\chapter*{\abstractname}
\bigskip
\begin{center}
\@title{}\par
\bigskip
\authors{}\par
\end{center}
\bigskip
}
\def\abstractname{ABSTRACT} % capitalized - MAL 5/3/2002


\newcommand{\acknowledgements}{\chapter*{\acknowledgementsname}}
\def\acknowledgementsname{Acknowledgments}
\newcommand{\preface}{\chapter*{\prefacename}}
\def\prefacename{Preface}
\newcommand{\introduction}{\chapter*{\introductionname}}
\def\introductionname{Introduction}
\newcommand{\vita}{\chapter*{\vitaname}}
\def\vitaname{Vita}

% Taken from amsbook.cls and modified.
% The references must use normal size font,
% and be single spaced with double space
% between items.
\renewenvironment{thebibliography}[1]{%
  \@xp\chapter\@xp*\@xp{\bibname}%
%  \normalfont\footnotesize            % commented out (MAL-5/4/2001)
  \labelsep .5em\relax
  \renewcommand\theenumiv{\arabic{enumiv}}\let\p@enumiv\@empty
  \list{\@biblabel{\theenumiv}}{%
    \settowidth\itemsep\@mainsize      % MAL - 5/9/2001
%    \singlespace                       % MAL - 5/9/2001
    \begin{singlespace}
    \settowidth\labelwidth{\@biblabel{#1}}%
    \leftmargin\labelwidth \advance\leftmargin\labelsep
    \usecounter{enumiv}}%
  \sloppy \clubpenalty\@M \widowpenalty\clubpenalty
  \sfcode`\.=\@m
}{%
  \def\@noitemerr{\@latex@warning{Empty `thebibliography' environment}}%
%  \endsinglespace                      % MAL - 5/9/2001
    \end{singlespace}
  \endlist
}
\def\bibname{References}


\def\appendices{\par
  \c@chapter\z@ \c@section\z@
  \let\chaptername\appendixname
  \def\thechapter{\@Alph\c@chapter}}

\def\appendix{\par
  \c@chapter\z@ \c@section\z@
  \let\chaptername\appendixname
  \def\thechapter{{}}
  \renewcommand\thesection{\@arabic\c@section}}


\doublespacing


\endinput
%%
%% End of file `nuthesis.cls'.

