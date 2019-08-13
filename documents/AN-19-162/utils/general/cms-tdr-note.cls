\def\fileversion{1.0.0}\def\filedate{21 December 2006}

\typeout{=======================================================}
\typeout{ }
\typeout{        CMS TDR/NOTE style. Version \fileversion}
\typeout{ }
\typeout{        Lucas Taylor,  \filedate}
\typeout{ }
\typeout{=======================================================}

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cms-tdr-note}[2006/12/21]
\usepackage{ifthen}%           Conditionals     

% switches to handle document styles: one of these must be specified
\newboolean{cms@tdr}
\setboolean{cms@tdr}{false} 
\newboolean{cms@note}
\setboolean{cms@note}{false}
\newboolean{cms@an}
\setboolean{cms@an}{false}
\newboolean{cms@in}
\setboolean{cms@in}{false}
\newboolean{cms@cr}
\setboolean{cms@cr}{false}

\DeclareOption{tdr}{\setboolean{cms@tdr}{true}\def\cms@style{tdr}} % default anyway
\DeclareOption{note}{\setboolean{cms@note}{true}\def\cms@style{note}}
\DeclareOption{an}{\setboolean{cms@an}{true}\def\cms@style{an}}
\DeclareOption{cr}{\setboolean{cms@cr}{true}\def\cms@style{cr}}
\DeclareOption{in}{\setboolean{cms@in}{true}\def\cms@style{in}}
% unknown options are passed on to report class or article class,
% depending on the cms@tdr switch (its value is known as its option is declared first)
\DeclareOption*{%
\ifthenelse{\boolean{cms@tdr}}%
{\typeout{Using report \CurrentOption\ option}\PassOptionsToClass{\CurrentOption,openright}{report}}%
{\typeout{Using article \CurrentOption\ option}\PassOptionsToClass{\CurrentOption,titlepage}{article}}%
{}}

\ProcessOptions\relax

\ifthenelse{\boolean{cms@note}}{\typeout{Using cms note sub-class}}{\typeout{Not using cms note sub-class}}{}
\ifthenelse{\boolean{cms@tdr}}
{\LoadClass{report}}
{\LoadClass{article}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  General style files 
%  -------------------
%
\usepackage{fancybox}%         e.g. for calling out requirements, etc.             
\usepackage{multicol}%         multiple columns (e.g. glossary)
%\usepackage{lineno}%           Does line numbers in margins (see \linenumbers command)

%\usepackage{cite}% 
\providecommand{\citen}[1]{\cite{#1}}%   

\usepackage{fancyhdr}%        modern version of fancyheadings
\usepackage{subfigure}%
\usepackage{xspace}% 

\usepackage{amsmath,amssymb}% for advanced math typesetting
\DeclareMathOperator{\sgn}{sgn}% Declare the sgn operator
\DeclareMathOperator{\sign}{sign}% for those who insist on sign instead

\usepackage{topcapt}% for captions at top, i.e., tables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  PDF and graphics settings 
%  -------------------------
%
\usepackage[pdftex]{graphicx,graphics}
\DeclareGraphicsExtensions{.pdf,.jpg,.jpeg}  % used by includegraphics. Note: .eps and .gif not supported by pdftex
\usepackage{lscape} % part of graphicx. used for rotating entire pages.
\usepackage{rotating}%  LT: added 3 Apr 2006 for Vincenzo 

% Use following line for electronic copy (coloured links etc)
%\usepackage[pagebackref=true,bookmarksnumbered,bookmarksopen,bookmarksopenlevel=1,colorlinks,linkcolor=magenta,citecolor=blue,urlcolor=red,plainpages=false,pdfpagelabels]{hyperref}
%\usepackage[bookmarksnumbered,bookmarksopen,bookmarksopenlevel=1,colorlinks,linkcolor=magenta,citecolor=blue,urlcolor=red,plainpages=false,pdfpagelabels]{hyperref}

% Use following line for paper (no coloured links etc)
\usepackage[bookmarksnumbered,bookmarksopen,bookmarksopenlevel=1,colorlinks=false,pdfborder=0,plainpages=false,pdfpagelabels]{hyperref}

\usepackage[numbers,sort&compress]{natbib}
\usepackage{hypernat}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
% Override the (rather clunky) default fonts of LaTex throughout
% --------------------------------------------------------------
%
\usepackage{palatino}%    Choose default roman font.  Others are times, pslatex, newcent, bookman, chancery
% can use mathpazo if no one uses \mathrm{\gamma}, etc., since it doesn't supply Greek in rm. Google "mathpazo mathrm greek"
\ifthenelse{\boolean{cms@tdr}}{}{\usepackage{mathpazo}}%    Matching math fonts (see http://www.math.uiuc.edu/~hartke/computer/latex/survey/survey.html)
\usepackage{helvet}%      Choose default sans serif

\providecommand{\chaptername}{Chapter} % not present in article class
\usepackage{sectsty}%     Change style of headings (should precede section redefinitions below)
\allsectionsfont{\sffamily}% Set sans serif for all headings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define the default parskip so it can be reset back to this after
% it is locally changed in other environments
%
\newlength{\parskipsaved}
\setlength{\parskipsaved}{0.5\baselineskip} 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fine tune definitions of sections   (from report.sty. See p 27, Latex companion)
% ---------------------------------
%
%------------------------start of chapter redefinition
\def\@makechapterhead#1{%
  \vspace*{30\p@}%  Was 50
  {\parindent \z@ \raggedright \normalfont
    \ifnum \c@secnumdepth >\m@ne
        \Large\sffamily\bfseries \@chapapp\space \thechapter % was huge
        \par\nobreak
        \vskip 10\p@ % was 20
    \fi
    \interlinepenalty\@M
    \huge\sffamily\bfseries #1\par\nobreak % Was Huge
    \vskip 40\p@ % was 40
  }}
\def\@schapter#1{\if@twocolumn
                   \@topnewpage[\@makeschapterhead{#1}]%
                 \else
                   \@makeschapterhead{#1}%
                   \@afterheading
                 \fi}
\def\@makeschapterhead#1{%
  \vspace*{30\p@}%
  {\parindent \z@ \raggedright
    \normalfont
    \interlinepenalty\@M
    \huge \sffamily \bfseries  #1\par\nobreak
    \vskip 40\p@
  }}
%------------------------end of chapter redefinition

\renewcommand\section{\@startsection 
{section}%
{1}%
{\z@}%
{-2.5ex \@plus -1ex \@minus -.2ex}%       Was: {-3.5ex \@plus -1ex \@minus -.2ex}%
{0.3ex \@plus.2ex}%                       Was: {2.3ex \@plus.2ex}
{\sffamily\Large\bfseries}}%              Was: \normalfont ...

\renewcommand\subsection{\@startsection
{subsection}%
{2}%
{\z@}%
{-1.0ex\@plus -1ex \@minus -.2ex}%       Was: {-3.25ex\@plus -1ex \@minus -.2ex}%  
{0.05ex \@plus .2ex}%                     Was:   {1.5ex \@plus .2ex}%
{\sffamily\large\bfseries}}%           Was: \normalfont ...

\renewcommand\subsubsection{\@startsection
{subsubsection}%
{3}%
{\z@}%
{-1.0ex\@plus -1ex \@minus -.2ex}%       Was: {-3.25ex\@plus -1ex \@minus -.2ex}
{0.05ex \@plus .2ex}%                     Was:   {1.5ex \@plus .2ex}%
{\sffamily\normalsize\bfseries}}%           Was: \normalfont ...

\renewcommand\paragraph{\@startsection
{paragraph}%
{4}%
{\z@}%
{-1.0ex\@plus -1ex \@minus -.2ex}%       Was: {3.25ex \@plus1ex \@minus.2ex}%             Was: 
{-2em}%                             Was:  If negative, means in-line. Space is then the indentation.
{\sffamily\normalsize\bfseries}}%           Was: \normalfont ...

\renewcommand\subparagraph{\@startsection
{subparagraph}%
{5}%
{\parindent}%
{-1.0ex\@plus -1ex \@minus -.2ex}%
{-2em}%
{\sffamily\normalsize\bfseries}}%           Was: \normalfont ...



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Formatting of table of contents, etc.
% -------------------------------------
%
\renewcommand\tableofcontents{%
\setlength{\parskip}{0.15\baselineskip}%%% LT: make lines tight 
  \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \chapter*{\contentsname
        \@mkboth{%
           \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
    \if@restonecol\twocolumn\fi
\setlength{\parskip}{\parskipsaved}%%% LT: reset to default value
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define page layout
% ------------------
%
\setlength{\textwidth}{156mm}
\setlength{\textheight}{230mm}
\setlength{\columnsep}{10mm}
\setlength{\columnseprule}{0.3pt}% if zero then there is no rule between 2 columns

\if@twoside
  \setlength{\oddsidemargin}{0mm}
  \setlength{\evensidemargin}{0mm}
\else
  \setlength{\oddsidemargin}{0mm}
  \setlength{\evensidemargin}{0mm}
\fi

\setlength{\footskip}{36pt}
\setlength{\headheight}{25pt}
\setlength{\headsep}{20pt}%     Distance from bottom of header to main body
\setlength{\topmargin}{0mm}

%
% Margin notes
% ------------
\setlength{\marginparpush}{\baselineskip}% space between successive notes   
\setlength{\marginparsep}{3mm}%            distance from main text
\setlength{\marginparwidth}{20mm}%         width of note

%
% Paragraph spacing / indentation
% ------------------------------- 
\setlength{\parindent}{0pt}
\setlength{\parskip}{\parskipsaved} 


\raggedbottom



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Float handling parameters for better placing of figures 
% -------------------------------------------------------
%
\renewcommand{\topfraction}{1.0}
\renewcommand{\dbltopfraction}{1.0}
\renewcommand{\bottomfraction}{1.0}
\renewcommand{\textfraction}{0.2}

\setlength{\intextsep}{10pt plus 2pt minus 4pt}       % Space between in-text figures and text.
\setlength{\floatsep}{10pt plus 2pt minus 4pt}        % Space between adjacent floats moved to top or bottom of text page.
\setlength{\textfloatsep}{10pt plus 2pt minus 4pt}    % Space between main text and floats at top or bottom of page.
\setlength{\dblfloatsep}{10pt plus 2pt minus 4pt}     % Same for double column
\setlength{\dbltextfloatsep}{10pt plus 2pt minus 4pt}

\setlength\abovecaptionskip{4pt plus 2pt minus 2pt}   % default is for caption at bottom. 
                                                      % use topcapt.sty: \topcaption for cap at top
\setlength\belowcaptionskip{0pt}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Miscellaneous style settings
% ----------------------------
%
% \bibliographystyle{unsrt}
\bibliographystyle{auto_generated}  % File generated by tdr script

\setcounter{secnumdepth}{4}   
\setcounter{tocdepth}{2}      

\renewcommand{\baselinestretch}{0.99}%   LT 2005-10-07   was 1.0
\normalsize% required to activate \baselinestretch

% \linenumbers   % shows line numbers in margins
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make \mychapter, \mysection, etc.
% These used to be more ``fancy'' but now revert to (pretty much)
% their default definitions.  LT, 12-oct-2005

\newcommand{\mychapter}[1]{\clearpage\chapter{#1}}
\newcommand{\mysection}[1]{\bigbreak\section{#1}}
\newcommand{\mysubsection}[1]{\bigbreak\subsection{#1}}
\newcommand{\mysubsubsection}[1]{\bigbreak\subsubsection{#1}}
\newcommand{\myparagraph}[1]{\bigbreak\paragraph{#1}}%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define page header and footer content
% -------------------------------------
%
% GOA 2006-12-26: use non-deprecated versions of fancyhdr macros 
\pagestyle{fancy}
\fancyhead[LE,RO]{\sffamily\bfseries{}\thepage}
\fancyhead[LO,RE]{\sffamily\small\rightmark}%  LT 2005-10-07 add \small's 

\fancyfoot{}% clear

% Redefine format of text in header lines (see. fancyhdr.sty)
% ---------------------------------------
% \def\MakeUppercase{}%       By default is:  \def\MakeUppercase{\uppercase}
% \def\MakeUppercase{\sc}%    By default is:  \def\MakeUppercase{\uppercase}
\def\MakeUppercase{\textbf}%     By default is:  \def\MakeUppercase{\uppercase}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define standard section headings
% --------------------------------
%  
\renewcommand{\contentsname}{Contents}
\renewcommand{\listfigurename}{Figures}
\renewcommand{\listtablename}{Tables}
\renewcommand{\bibname}{References}
\renewcommand{\indexname}{Index}
\renewcommand{\chaptername}{Chapter}
\renewcommand{\appendixname}{Appendix}
\renewcommand{\abstractname}{Abstract}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Re-New itemize environment (see e.g. P61 of LaTeX companion)
% --------------------------
%

\renewenvironment{itemize}%
{\begin{list}%
{$\bullet$}%
{\setlength{\topsep}{0.0ex}%        
 \setlength{\partopsep}{0.0ex}%        
 \setlength{\parsep}{0.5\parskip}% 
 \setlength{\itemsep}{0.0ex}% 
 \setlength{\leftmargin}{5ex}%  
 \setlength{\rightmargin}{5ex}%   
}}%
{%                                
\end{list}%
}%

% itemize same as itemize: for backwards compatibility
\newenvironment{myitemize}%
{\begin{list}%
{$\bullet$}%
{\setlength{\topsep}{0.0ex}%        
 \setlength{\partopsep}{0.0ex}%        
 \setlength{\parsep}{0.5\parskip}% 
 \setlength{\itemsep}{0.0ex}% 
 \setlength{\leftmargin}{5ex}%  
 \setlength{\rightmargin}{5ex}%   
}}%
{%                                
\end{list}%
}%

\newcommand{\marginlabel}[1]{\mbox{}\marginpar{\raggedleft\hspace*{0pt}#1}}%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\glossname{Glossary}
\providecommand{\GLOS}[2]{\makebox[0.110\textwidth][l]{\bfseries\mathversion{bold}#1}
\makebox[0.025\textwidth][l]{~}
\parbox[t]{0.32\textwidth}{#2\vspace*{-2.75mm}\\
}}
\newenvironment{theglossary}%
{\chapter*{\glossname}
\@mkboth{\MakeUppercase\glossname}{\MakeUppercase\glossname}%
% this is to get the pdfbookmark correct. should have come through @schapter but didn't?
  \begingroup
    \let\@mkboth\@gobbletwo
    \Hy@GlobalStepCount\Hy@linkcounter
    \xdef\@currentHref{\glossname*.\the\Hy@linkcounter}%
    \Hy@raisedlink{\hyper@anchorstart{\@currentHref}\hyper@anchorend}%
  \endgroup
\addcontentsline{toc}{chapter}{\glossname}%
\setlength{\columnseprule}{0.2mm}
\setlength{\multicolsep}{7mm}
%
\vspace*{-5mm}
\begin{multicols}{2}
\small%
\raggedcolumns%
\noindent}
{%
\end{multicols}
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Colour Plates environment
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\newcommand{\colourplates}{CP \@arabic\c@figure}
\def\CMScpname{Colour Plates}
\newenvironment{CMScolourplates}%
{\chapter*{\CMScpname}
\@mkboth{\MakeUppercase\CMScpname}{\MakeUppercase\CMScpname}%
% this is to get the pdfbookmark correct. should have come through @schapter but didn't?
  \begingroup
    \let\@mkboth\@gobbletwo
    \Hy@GlobalStepCount\Hy@linkcounter
    \xdef\@currentHref{CP*.\the\Hy@linkcounter}%
    \Hy@raisedlink{\hyper@anchorstart{\@currentHref}\hyper@anchorend}%
  \endgroup
\addcontentsline{toc}{chapter}{\CMScpname}%
\cleardoublepage
\begingroup
\renewcommand\thefigure{\colourplates}
}
{\endgroup}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modified bibliography processing (from article class)
%
\def\refname{References}
\renewenvironment{thebibliography}[1]
     {\chapter*{\refname}%
      \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}%
%
%      Note:\\ CHEP04 documents on InDiCo are available at \\%
%     \indent\url{http://indico.cern.ch/confAuthorIndex.py?confId=0}\\
%      CMS Internal Notes are available upon approved request from the CMS Secretariat.\\%
%     CMS Notes are available at \url{http://cmsdoc.cern.ch/docall.shtml} unless otherwise noted.
%
    \subsection*{Notes:}%
% this is to get the pdfbookmark correct. should have come through @schapter but didn't?
  \begingroup
    \let\@mkboth\@gobbletwo
    \Hy@GlobalStepCount\Hy@linkcounter
    \xdef\@currentHref{\refname*.\the\Hy@linkcounter}%
    \Hy@raisedlink{\hyper@anchorstart{\@currentHref}\hyper@anchorend}%
  \endgroup    
    \addcontentsline{toc}{chapter}{\refname}%     
    \begin{itemize}
%%DEA     \item[(a)] CMS Internal Notes and Analysis Notes are available upon approved request from the CMS Secretariat.\\%
     \item[(a)] CMS Notes are available at \url{http://cms.cern.ch/iCMS/} unless otherwise noted.\\%
     \item[(b)]  References marked \textbf{doi} should be prefixed with \url{http://dx.doi.org/}.%
    \end{itemize}
\rule{\linewidth}{1.0pt}
    %
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Expand LOF width so we can deal with Figure 10.10; this increases the label for tables as well
\renewcommand{\l@figure}{\@dottedtocline{1}{1.5em}{3.3em}}
\renewcommand{\l@section}{\@dottedtocline{1}{1.5em}{2.5em}}  %last arg was 2.3em
\renewcommand{\l@subsection}{\@dottedtocline{2}{3.8em}{3.5em}}  %last arg was 3.2em

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   NOTES code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ifthenelse{\boolean{cms@tdr}}{
}
{
\typeout{Defining CMS Note commands}
\pagestyle{empty}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% title page %%%%%%%%%%%%%%%%%%%%%%%%%
% author/address from class 'aps2005' for AP-S/URSI 2005 conference 
\newcounter{author}
\renewcommand{\author}[2][]{
   \stepcounter{author}
   \@namedef{author@\theauthor}{#2}
   \@namedef{authorlabel@\theauthor}{#1}
}

\newcounter{address}
\newcommand{\address}[2][]{
   \stepcounter{address}
   \@namedef{address@\theaddress}{#2}
   \@namedef{addresslabel@\theaddress}{#1}
}

\newcommand{\alsep}{and}

% uppercase the cms@style as note type: see TeXbook exercise 20.19 for uppercase tecnique
\newcommand{\cmsNoteHeader}[1]{\edef\next{\def\noexpand\cmsSTYLE{\cms@style}}\uppercase\expandafter{\next}\def\cmsNUMBER{#1}}
% Abstract: now a command, not an environment.
\newsavebox{\CMSabstractbox}
\savebox{\CMSabstractbox}{Please replace the default abstract using the abstract command.}
\renewcommand{\abstract}[1]{\savebox{\CMSabstractbox}{\parbox[t]{0.8\textwidth}{#1}}}
\newenvironment{CMSabstract}{\@beginparpenalty\@lowpenalty
    \begin{center}%
    \vskip1.5em\bfseries \abstractname
    \@endparpenalty\@M
    \end{center}
    \setlength{\parindent}{0pt}%
    \setlength{\parskip}{0.5\baselineskip}%
    \ignorespaces}
{\par}
% Definition of maketitle mostly from article.cls
\renewcommand\maketitle{\begin{titlepage}%
  \pagestyle{empty}
  \setcounter{page}{1}\pagenumbering{roman}
  \let\footnotesize\small
  \let\footnoterule\relax
  \let \footnote \thanks
  \topmargin -3.5cm
  \parskip 0pt
Available on the \href{http://cms.cern.ch/iCMS/jsp/iCMS.jsp?mode=single&part=publications}{CMS information server}
\hfill{\Large\textbf{CMS \cmsSTYLE\ \cmsNUMBER}}
\vskip -2em
\begin{center}\includegraphics{cms_\cmsSTYLE}\end{center}
\vskip -0.5em%
{\hspace\fill \large \@date \par}%       % Set date in \large size.
  \begin{center}%
    {\LARGE \@title \par}%
    \vskip 3em%
%    {\large
     \lineskip .75em%
     %%%%%%%%%% author block %%%%%%%%%%
    \newcounter{c@author}
    \newcounter{c@tmp}
    \ifthenelse{\value{author}=2}{%
      \newcommand{\liand}{ and }}{%
      \newcommand{\liand}{, and }}
    % Next add the authors
    \ifthenelse{\value{address}<2}{%
      %No address references required if there is only one
      \@nameuse{author@1}%
      \stepcounter{c@author}%
      \whiledo{\value{c@author}<\value{author}}{%
        \setcounter{c@tmp}{\value{author}}%
        \addtocounter{c@tmp}{-\value{c@author}}%
        \ifthenelse{\value{c@tmp}=1}{%
          \renewcommand{\alsep}{\liand}}{\renewcommand{\alsep}{, }}%
        \stepcounter{c@author}\alsep \@nameuse{author@\thec@author}}\\%
    }
    {%Add address references after the author's name
      \@nameuse{author@1}${}^{\ref{\@nameuse{authorlabel@1}}}$%
      \stepcounter{c@author}%
      \whiledo{\value{c@author}<\value{author}}{%
      \setcounter{c@tmp}{\value{author}}%
      \addtocounter{c@tmp}{-\value{c@author}}%
      \ifthenelse{\value{c@tmp}=1}{%
        \renewcommand{\alsep}{\liand}}{\renewcommand{\alsep}{, }}%
      \stepcounter{c@author}\alsep \@nameuse{author@\thec@author}%
        ${}^{\ref{\@nameuse{authorlabel@\thec@author}}}$%
      }
    }
    \vspace{0.2ex}

    \ifthenelse{\value{address}>0}{%
      \ifthenelse{\value{address}=1}{
        %Output the address without any reference symbols
       %%        {\slshape\small \@nameuse{address@1}}
        {\@nameuse{address@1}}
      }
      {%Output the addresses as an enumerated list
        \newcounter{c@address}

        \begin{center}
        \whiledo{\value{c@address}<\value{address}}
        {
          \refstepcounter{c@address}
            ${}^\thec@address$\,%
              \label{\@nameuse{addresslabel@\thec@address}}%
              \@nameuse{address@\thec@address}\\ %
        }
        \end{center}
      
      } % end of the address creation ifthenelse block
    }
    {
      \relax
    }
  \end{center}
  %%%%%%%%% end author block %%%%%%%%%%
%  }
  \par
  \begin{CMSabstract}
  \centering
  \usebox{\CMSabstractbox}
  \end{CMSabstract}
  %\@thanks
  \vfil\null
  \end{titlepage}%
  \cleardoublepage
  \pagestyle{fancy}
  \setcounter{page}{1}\pagenumbering{arabic}
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
  % from general tdr layout
\setlength{\footskip}{36pt}
\setlength{\headheight}{25pt}
\setlength{\headsep}{20pt}%     Distance from bottom of header to main body
\setlength{\topmargin}{0mm}
}
\renewenvironment{thebibliography}[1]
     {\section*{\refname}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}{}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Miscellaneous macros for diagnostics and options
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% macro to set bw or colour figures: \combinedfigure{switches for includegraphics}{bw figure name}{colour figure name}
\newboolean{bwfigurespreferredswitch}
\setboolean{bwfigurespreferredswitch}{false}
%\setboolean{bwfigurespreferredswitch}{true}
\providecommand{\combinedfigure}[3]{%
\ifthenelse{\boolean{bwfigurespreferredswitch}}%
{\includegraphics[#1]{#2}}%
{\includegraphics[#1]{#3}}%
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  define an annotation command for pdf: most from CGvanderLaan at hetnet.nl
\providecommand\annotation[3]%Note title, note text
{\pdfannot width 1cm height 0.5cm depth 0cm
  {/Subtype/Text
  /Open false                 %Open by default in viewer?
  /T (#1)            %Title, shows up as Reviewer Name in comment filter
  /C [#2]                %Colour
  /Contents (#3)             %Annotation text
  }
}
%  and use it to call out included files: use 'perl -i -p -e 's/input\s*{/tracinginput{/g' */*.tex' first (*/*/*.tex, etc to recurse)
\def\tracinginputswitch{false}
\providecommand{\tracinginput}[1]{\ifthenelse{\equal{\tracinginputswitch}{true}}
{\annotation{Input File}{0.9 0.9 0}{#1}\input{#1}}
{\input{#1}}
}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\def\tracinggraphicsswitch{false}
\makeatletter\let\OldGinclude@graphics\Ginclude@graphics\renewcommand{\Ginclude@graphics}[1]{\ifthenelse{\equal{\tracinggraphicsswitch}{true}}{\annotation{Graphics File}{0 0 .9}{#1}\typeout{<789FIG  #1 \thefigure>}\OldGinclude@graphics{#1}}{\OldGinclude@graphics{#1}}}\makeatother
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newboolean{colourplateswitch}
\setboolean{colourplateswitch}{true}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\endinput
