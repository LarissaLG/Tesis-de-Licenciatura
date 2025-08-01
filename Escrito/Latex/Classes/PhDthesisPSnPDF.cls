%-------------------------------------------------------------------
% Modification for the Science Faculty at UANM by Jonathan Urrutia
% Modificaciones para FI-UNAM por Jesús Velázquez y Marco Ruiz, 2014
% PhDthesis version 2 style file by Jakob Suckale, 2007
% based on CUEDthesis version 1 by Harish Bhanderi 2002
%-------------------------------------------------------------------

%-------------------------- identification ---------------------
%\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{Latex/Classes/PhDthesisPSnPDF}[2007/09/06 v2 PhD thesis class]


%If you want to use a Book style document then uncomment the following 3 lines and comment the above 8 report style lines:
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ProcessOptions\relax
\ifx\pdfoutput\undefined
   \LoadClass[dvips, letterpaper]{book}  % 
\else
   \LoadClass[pdftex, letterpaper]{book}
\fi

%-------------------------------------------------------------------------------
%                                   PACKAGES 		                                    |
%-------------------------------------------------------------------------------

\usepackage{lmodern}                                                        %Latin Modern font
%\usepackage{pslatex}						%Times New Roman font 
\usepackage{multicol}						% for pages with multiple text columns, e.g. References
	\setlength{\columnsep}{20pt} 				% space between columns; default 10pt quite narrow
\usepackage{multirow}
\usepackage{listings}                    				% Allows to add coding lines 
\usepackage{emptypage}						% Eliminates headers and footnotes in empty pages
\usepackage{makeidx}						% Add glosaries
%\usepackage[style=list,toc,number=none]{glossary}
\raggedbottom								% Makes LaTeX send tha white spaces to the background
\usepackage[bottom]{footmisc}			% for the footnotes; perpage -> inicia la numeración en la página misma
\usepackage{fancyhdr}						% for better header layout
\usepackage{eucal}
\usepackage{ifthen}
\usepackage[nottoc]{tocbibind} 					% correct page numbers for bib in TOC, nottoc suppresses an entry for TOC itself
%\usepackage{nextpage}
\usepackage{titlesec}
\usepackage{imakeidx}%for indexes
\makeindex[intoc]
\usepackage{tocloft}
\usepackage{xcolor}
\definecolor{myTurquoise}{HTML}{12B3C5} 

\setlength\cftparskip{.6em}
\setlength\cftbeforechapskip{0pt}
	%---------------------- Figurbackref change legend latex bibtexes ----------------------------
		\usepackage[margin=10pt,font=scriptsize,labelfont=bf, ]{caption}  
		\usepackage[labelformat=simple]{subcaption} % labelformat=simple   O brace
		\usepackage{caption}
		\usepackage{floatrow}

		\floatsetup[figure]{style=plain,subcapbesideposition=top}

%			\captionsetup[subfigure]{labelformat=simple}
		%\usepackage[bf,SL,BF]{subfigure}         
                     % \usepackage{float}
		%	\floatsetup[subfigure]{style=plain,heightadjust=object,
		%					capbesideposition={left,top},
		%					capbesidesep=space}
			%\DeclareCaptionSubType[alph]{figure}
			%\captionsetup[subfigure]{labelformat=brace,justification=centerlast}
		\usepackage[export]{adjustbox}
		\usepackage{wrapfig}                    			% to include figure with text wrapping around it
		\usepackage{subfig}
		\captionsetup[subfigure]{labelformat=simple, labelsep=none}
		\renewcommand{\thesubfigure}{\alph{subfigure})}

	%--------------------- Color/equations -related stuff
		\usepackage[dvipsnames]{xcolor}
			%UNAM color palette
			\definecolor{UNAMblue}{RGB}{0,60,113}		% Blue pantone  541 (C)
			\definecolor{UNAMgold}{RGB}{234,221,150}		% Gold pantone  460  (C) 
			\newcommand\myshade{85}
			\colorlet{mylinkcolor}{violet}
			\colorlet{mycitecolor}{YellowOrange}
			\colorlet{myurlcolor}{Aquamarine}
		\usepackage{empheq}
		\usepackage[most]{tcolorbox}
			\tcbset{colframe=myTurquoise!90!blue,  % color base para el marco
    colback=myTurquoise!20!white,    colbacktitle = black, boxrule=1.5pt,  
			fonttitle=\bfseries , standard jigsaw,  opacityback=0.95, opacitybacktitle = .5,center}% sharp corners = downhill, boxrule = .75mm }

\newcommand{\coloredeq}[2]{\begin{empheq}[box= \fcolorbox{myTurquoise!90!blue}{myTurquoise!20!white},{\fboxsep=15pt\fbox}]{align}\label{#1}#2\end{empheq}}
				
	%----------------------- Language and references ------------------------------
		\usepackage[spanish,mexico]{babel}			% Comment for english
			\addto\captionsspanish{\renewcommand{\figurename}{Fig.}}		
		\usepackage[T1]{fontenc}			 		% 8-bits fonts
		\usepackage[utf8]{inputenc}           			% Complete latin alfabet: Accents and so on
		\spanishdecimal{.}
		\usepackage{microtype}
		\usepackage[bibstyle = trad-abbrv, 
					citestyle = numeric-comp,
					backend = bibtex,
					sorting = none,
					backref=true]{biblatex} %style = trad-abbrv
		\DeclareFieldFormat[article]{volume}{\mkbibbold{#1}}
		\AtEveryBibitem{\clearfield{urldate}}
		\AtEveryBibitem{\clearfield{url}}
		\AtEveryBibitem{\clearfield{isbn}}
		\AtEveryBibitem{\clearfield{month}}
		\AtEveryBibitem{\clearfield{day}}
		\AtEveryBibitem{\clearfield{issn}}

		\DefineBibliographyStrings{spanish}{%
			  backrefpage = {citado en la pág.},% originally "cited on page"
			  backrefpages = {citado en las págs.},% originally "cited on pages"
		}	
	
		\addbibresource{6-References/references.bib}
		\addto\captionsspanish{\renewcommand{\figurename}{Fig.}}

	%------------------------- Maths --------------------------------------
		\usepackage{amssymb, amsmath, amsbsy, amsfonts}	% Math symbols and so on
		\usepackage{mathrsfs}						
		\usepackage{mathdots}                    			%  \iddots
		\usepackage{mathtools}
		\usepackage{physics}


%---------------------------------------------------------------------------------
%                                   CHAPTERS FORMAT 		                       |
%---------------------------------------------------------------------------------

		\newcommand{\bigrule}{\titlerule[0.5mm]}

\titleformat{\chapter}[display] 		% Changes chapter forms
{\Huge} 								% Sets the font size (and format if wished e.g. {\Huge\bfseries}
{  
\titlerule 								% Horizontal line
\filleft 								% Fills the left size od the line so the text becomes right aligned
\Large\chaptertitlename\ 				% Changes to \Large insted of \Huge (but why?)
\Large\thechapter}						% Changes to \Large insted of \Huge (but why?)
{0mm} 									% Minimum size between textbody and label (??)
{\filleft} 								% Right aligning
[\vspace{0.5mm} \bigrule] 				% Bold horizontal line beneath the chapter name


%:-------------------------- Glossary/Abbrev./Symbols ----------------------- 
%\usepackage[intoc]{nomencl} % load nomencl extension; include in TOC
%\nomrefpage % to include page numbers after abbrevations
%\renewcommand{\nomname}{Glossary} % rename nomenclature
%\renewcommand{\nomlabel}[1]{\textbf{#1}} % make abbreviations bold
%\makenomenclature % used to be \makeglossary
%\newcommand{\g}{\footnote{For all abbreviations see the glossary on page \pageref{nom}.}} % type "\g" to refer to glossary

% used to be for sorting into categories:
%\renewcommand\nomgroup[1]{%
%  \ifthenelse{\equal{#1}{A}}{%
%   \item[\textbf{Roman Symbols}] }{%             A - Roman
%    \ifthenelse{\equal{#1}{G}}{%
%     \item[\textbf{Greek Symbols}]}{%             G - Greek
%      \ifthenelse{\equal{#1}{R}}{%
%        \item[\textbf{Superscripts}]}{%              R - Superscripts
%          \ifthenelse{\equal{#1}{S}}{%
%           \item[\textbf{Subscripts}]}{{%             S - Subscripts
%	    \ifthenelse{\equal{#1}{X}}{%
%	     \item[\textbf{Other Symbols}]}{{%    X - Other Symbols
%	    \ifthenelse{\equal{#1}{Z}}{%
%	     \item[\textbf{Acronyms}]}%              Z - Acronyms
%              			{{}}}}}}}}}}


%---------------------------------------------------------------------------------
%                                  PDF/PS set up 				hyperref package		| 
%---------------------------------------------------------------------------------

%if you use a macTeX 2008 or later, use the ifpdf package
\usepackage{ifpdf}

%if you use an older version, uncomment these lines:
%\newif \ifpdf
%    \ifx \pdfoutput \undefined
        % for running latex
%        \pdffalse
%    \else
        % for running pdflatex
%        \pdfoutput = 1    % positive value for a PDF output
                          % otherwise a DVI
%        \pdftrue
%\fi

\ifpdf

    \usepackage[ pdftex, plainpages = false, pdfpagelabels, 
                 pdfpagelayout = OneColumn, 
                 bookmarks,
                 bookmarksopen = true,
                 bookmarksnumbered = true,
                 breaklinks = true,
                 linktocpage,
                 colorlinks = true,
                 linkcolor =  blue,
                 urlcolor  = magenta,
                 citecolor = red,
                 anchorcolor = green,
				hyperfootnotes = false,
                 hyperindex = true,
                 hyperfigures
                 ]{hyperref} 
%	\hypersetup{ hidelinks = true, }% Uncomment for a print version, where the hyperlinks have no colors.
    \usepackage{graphicx}   
    \DeclareGraphicsExtensions{.png, .jpg, .pdf}

    \pdfcompresslevel=9
    %\graphicspath{{ThesisFigs/PNG/}{ThesisFigs/PDF/}{ThesisFigs/}}
\else
    \usepackage[ dvips, 
                 bookmarks,
                 bookmarksopen = true,
                 bookmarksnumbered = true,
                 breaklinks = true,
                 linktocpage,
%                 pagebackref,
                 colorlinks = true,
                 linkcolor = blue,
                 urlcolor  = magenta,%blue,
                 citecolor = green,%lgray,
                 anchorcolor = dgreen,
                 hyperindex = true,
                 hyperfigures
                 ]{hyperref}

    %\usepackage{epsfig}
    \usepackage{graphicx}
    \DeclareGraphicsExtensions{.eps, .ps}
    %\graphicspath{{ThesisFigs/EPS/}{ThesisFigs/}}
\fi

%----------------------------------------------------------------------------------
%                               Page layout							|  NO MODIFICATIONS
%----------------------------------------------------------------------------------

\ifpdf				%A4 settings
   \pdfpageheight=297mm
   \pdfpagewidth=210mm
\else
   \setlength{\paperheight}{297mm}
   \setlength{\paperwidth}{210mm}
\fi

\setlength{\hoffset}{0.00cm}
\setlength{\voffset}{0.00cm}

%%: Uncomment this secion for two-sided printing
%% ------------------------------
%\setlength{\oddsidemargin}{1cm}  
%\setlength{\evensidemargin}{-.5cm} 
%\setlength{\topmargin}{-.75cm}  
%\setlength{\headheight}{0.cm}		
%\setlength{\headsep}{1cm}			
%\setlength{\textheight}{22.5cm}   
%\setlength{\textwidth}{16cm}  		
%\setlength{\marginparsep}{1mm}
%\setlength{\marginparwidth}{3cm}		
%%\setlength{\footskip}{2.36cm}
%\setlength{\footskip}{0.7cm}
%\setlength{\parskip}{1 em}	
%\setlength{\parindent}{3em}   

%%: Uncomment this secion for one-sided printing or digital version
%% ------------------------------
\setlength{\oddsidemargin}{.5cm}  
\setlength{\evensidemargin}{.5cm}
\setlength{\topmargin}{-.75cm}  
\setlength{\headheight}{0.cm}		
\setlength{\headsep}{1cm}			
\setlength{\textheight}{22.5cm}   
\setlength{\textwidth}{16cm}  		
\setlength{\marginparsep}{1mm}
\setlength{\marginparwidth}{3cm}		
%\setlength{\footskip}{2.36cm}
\setlength{\footskip}{0.7cm}
\setlength{\parskip}{1 em}	
\setlength{\parindent}{3em}   

	%------------------------- Fancy page layout options --------------------------------------

		\pagestyle{fancy}	
		\renewcommand{\chaptermark}[1]{\markboth{\MakeUppercase{\thechapter. #1 }}{}}
		\renewcommand{\sectionmark}[1]{\markright{\thesection\ #1}}

		\fancyhf{}
		\fancyhead[LO]{\rightmark}
		\fancyhead[LE]{\leftmark}
		\fancyfoot[LE, RO]{\thepage }
		\renewcommand{\headrulewidth}{0.5pt}
		\renewcommand{\footrulewidth}{0.1pt}

		\addtolength{\headheight}{0.5pt}
		\fancypagestyle{plain}{
			  \fancyhead{}
			  \renewcommand{\headrulewidth}{0pt} }
	
		\setlength{\skip\footins}{.5cm} 			% gap between text and footer




	%------------------------- Cover page layout options --------------------------------------

		% starts roman page numbering until chapter 1
		% important to avoid two pages numbered 1 and 2 which may cause bad links
		% bug: cover i + back side ii and then numbering restarts with i; should be iii
		\renewcommand{\thepage}{\roman{page}}

		%At the end of thissocument the arabic nubers are set

		%------------------------- Macro definitions --------------------------------------
		% These macros are used to declare arguments needed for the construction of the title page and other preamble.
		
			%------------------------- Personal info --------------------------------------
		
				\newcommand{\@degreedate}{}     			% The year and term the degree will be officially conferred		
				\newcommand{\degreedate}[1]{\renewcommand{\@degreedate}{#1}}
	
				\newcommand{\@degree}{}				% The full (unabbreviated) name of the degree
				\newcommand{\degree}[1]{\renewcommand{\@degree}{#1}}
		
				\newcommand{\@dept}{} 					% The name of your department(eg. Engineering, Maths, Physics)
				\newcommand{\dept}[1]{\renewcommand{\@dept}{#1}}

				\newcommand{\@college}{}				% The name of your college (eg. King's)
				\newcommand{\college}[1]{\renewcommand{\@college}{#1}}

				\newcommand{\@university}{}				% The name of your University
				\newcommand{\university}[1]{\renewcommand{\@university}{#1}}

				\newcommand{\@escudouniversidad}{}		% Defining the univertity crest
				\newcommand{\escudouniversidad}[1]{\renewcommand{\@escudouniversidad}{#1}}

				\newcommand{\@escudofacultad}{}			% Defininf the faculty/school crest
				\newcommand{\escudofacultad}[1]{\renewcommand{\@escudofacultad}{#1}}

			%------------------------- Asessors' info--------------------------------------

				\newcommand{\@director}{}
				\newcommand{\director}[1]{\renewcommand{\@director}{#1}}

				\newcommand{\@ctutoruno}{}
				\newcommand{\ctutoruno}[1]{\renewcommand{\@ctutoruno}{#1}}

				\newcommand{\@ctutordos}{}
				\newcommand{\ctutordos}[1]{\renewcommand{\@ctutordos}{#1}}

				\newcommand{\@presidente}{}
				\newcommand{\presidente}[1]{\renewcommand{\@presidente}{#1}}														

				\newcommand{\@vocal}{}
				\newcommand{\vocal}[1]{\renewcommand{\@vocal}{#1}}

				\newcommand{\@secretario}{}
				\newcommand{\secretario}[1]{\renewcommand{\@secretario}{#1}}
								
				\newcommand{\@supuno}{}
				\newcommand{\supuno}[1]{\renewcommand{\@supuno}{#1}}												

				\newcommand{\@supdos}{}
				\newcommand{\supdos}[1]{\renewcommand{\@supdos}{#1}}										
	
				\newcommand{\@institucion}{}
				\newcommand{\institucion}[1]{\renewcommand{\@institucion}{#1}}
				
				\newcommand{\@pags}{}
				\newcommand{\pags}[1]{\renewcommand{\@pags}{#1}}
				
				\newcommand{\@student}{}
				\newcommand{\student}[1]{\renewcommand{\@student}{#1}}							

				\newcommand{\@programa}{}
				\newcommand{\programa}[1]{\renewcommand{\@programa}{#1}}

				\newcommand{\@campo}{}
				\newcommand{\campo}[1]{\renewcommand{\@campo}{#1}}

			%------------------------- Place info --------------------------------------

				\newcommand{\@lugar}{}
				\newcommand{\lugar}[1]{\renewcommand{\@lugar}{#1}}

			%------------------------- Covers' definition --------------------------------------
				\newif\ifportada
				\newif\ifposgrado
				\newif\ifcomite

			%------------------------- Metadata and others--------------------------------------

				\newcommand{\submissiontext}{This dissertation is submitted for the degree of }  % Submission Text

				\newcommand{\@keywords}{}	% These keywords will appear in the PDF meta-information called `pdfkeywords`.
				\newcommand{\keywords}[1]{\renewcommand{\@keywords}{#1}}

				\newcommand{\@subject}{}
				\newcommand{\subject}[1]{\renewcommand{\@subject}{#1}}

			% These macros define an environment for front matter that is always  single column even in a double-column document.
 				\newenvironment{alwayssingle}{
			       \@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
			       \else\newpage\fi
				}
				{\if@restonecol\twocolumn\else\newpage\fi}

%---------------------------------------------------------------------------------
%                                  UNAM COVER 						| 
%---------------------------------------------------------------------------------

	%------------------- lines commands
		\newcommand{\colorvrule}[3]{\begingroup\color{#1}\vrule width#2 height#3\endgroup}
		\newcommand{\colorhrule}[2]{   \begingroup\color{#1}\hrule height#2		\endgroup}		 %color,  width

\renewcommand{\maketitle}{
\ifposgrado
	\begin{alwayssingle}
	\thispagestyle{empty}	
	\begin{center}
	\includegraphics[height=3cm]{Latex/Classes/Escudos/UNAM_crest_blue}\\[8mm]
	{\bfseries UNIVERSIDAD~NACIONAL~AUTÓNOMA~DE~MÉXICO}\\
	{\MakeUppercase {\@programa}}\\
	{\MakeUppercase {\@campo}}\\[20mm]
	{\Large \scshape {\@title}} \\[20mm]
	TESIS\\
	QUE PARA OPTAR POR EL GRADO DE:\\
	{\Large \scshape {\@degree}}\\[20mm]
	PRESENTA:\\
	{\Large \scshape {\@author}}\\[18mm]
	DIRECTOR DE TESIS\\
	{\Large \scshape {\@director}}\\
	\ifcomite
		COMITÉ TUTOR\\
		{\Large \scshape {\@ctutoruno}}\\
		{\Large \scshape {\@ctutordos}}\\[10mm]
		{\MakeUppercase {\@lugar}}{, }{\MakeUppercase {\@degreedate}}
	\else
		\vspace{20mm}
		{\MakeUppercase {\@lugar}}{, }{\MakeUppercase {\@degreedate}}
	\fi	
	\end{center}
	\end{alwayssingle}

	\newpage
	\thispagestyle{empty}~
	
		
	\newpage
	\begin{alwayssingle}
	\thispagestyle{empty}
	\vspace{30mm}

	\noindent JURADO ASIGNADO:
	\vspace{8mm}

	\begin{tabbing}
	Presidente: \ \ \ \ \ \ \ \= {\@presidente}\\ [3mm]
	Secretario: \> {\@secretario}\\ [3mm]
	Vocal: \> {\@vocal}\\ [3mm]
	1er. Suplente: \> {\@supuno}\\ [3mm]
	2o. Suplente: \> {\@supdos}\\ [12mm]
	\end{tabbing}
	La tesis se realizó en {\@institucion}.
	\vspace{15mm}

	\begin{center}
	DIRECTOR DE TESIS:\\[4mm]
	{\@director}\\[20mm]
	---------------------------------------------------------\\
	\end{center}
	%\newpage
	\vfill
	\end{alwayssingle}
	

\else			% Meaning bachelor's work
	\begin{alwayssingle}
	\begin{center}
	\thispagestyle{empty}

	%------------------------- Crests and left bar
		\hskip -2cm
		\begin{minipage}[c][10cm][s]{3cm}
		\begin{center}
			\ifportada
				\includegraphics[height=2.6cm]{Latex/Classes/Escudos/UNAM_crest_blue}\\[10pt]
				\hskip2pt
%				\colorvrule{UNAMgold}{1pt}{13cm}
%				\hskip -0.4mm
				\colorvrule{UNAMblue}{2.5pt}{13cm}
				\hskip -0.4mm
				\colorvrule{UNAMgold}{1pt}{13cm}\\[10pt]
				\includegraphics[height=2.6cm]{Latex/Classes/Escudos/FC_crest_blue}
			\else
				\includegraphics[height=2.6cm]{Latex/Classes/Escudos/UNAM_crest_black}\\[10pt]
			        \hskip2pt
				\vrule width2.5pt height13cm
				\hskip1mm
			        \vrule width1pt height13cm\\[10pt]
	       			\includegraphics[height=2.6cm]{Latex/Classes/Escudos/FC_crest_black}
			\fi
		\end{center}
		\end{minipage}\quad
	%------------------------- Upper bar and titles
		\begin{minipage}[c][9.5cm][s]{11cm}
		 \begin{center}
		%------------------------- Upper bar
			{\large \scshape Universidad Nacional Autónoma de México}
    			\vspace{.3cm}
			\ifportada
				%\hrule height2.5pt
%				\colorhrule{UNAMgold}{1pt}
%				\vspace{.1cm}				
				\colorhrule{UNAMblue}{2.5pt}
				\vspace{.1cm}
				%\hrule height1pt
				\colorhrule{UNAMgold}{1pt}	
			\else
				\hrule height2.5pt
				\vspace{.1cm}
				\hrule height1pt
			\fi
			\vspace{.3cm}
			{\scshape  FACULTAD DE CIENCIAS}
			\vspace{.3cm}
	
		%------------------------- Work title
			%\vspace{3cm}
			\vspace{2.5cm}	
			{\Large \scshape {\@title}}
			\vspace{3cm}

		%------------------------- Work type
			\makebox[8cm][s]{\Huge T E S I S}\\[8pt]
			QUE PARA OPTAR POR EL GRADO DE:\\[5pt]
			\textbf{{\@degree}}\\[40pt]
			PRESENTA:\\[5pt]
			\textbf{{\@author}}
			 \vspace{1.5cm}

		%-------------------------Tutor
			{\small DIRECTOR DE TESIS:\\ {\@director}}
			\vspace{0.5cm}
			%{\small CODIRECTOR DE TESIS:\\ en caso de que aplique}
			\vspace{1cm}

		%-------------------------Place
	    		{\@lugar}{, }{\@degreedate}
		
		\end{center}
		\end{minipage}

	\end{center}
	\null\vfill
	\end{alwayssingle}

	\newpage
	\thispagestyle{empty}~
	
	\newpage
	\begin{alwayssingle}
	\thispagestyle{empty}
	\vspace{30mm}

	\noindent 1. Datos del alumno\\ \@student \\[1em]
	\noindent 2. Datos del tutor \\ \@secretario \\[1em]
	\noindent 3. Datos del sinodal 1\\ \@presidente \\[1em]
	\noindent 4. Datos del sinodal 2 \\ \@vocal \\[1em]
	\noindent 5. Datos del sinodal 3\\ \@supuno \\[1em]
	\noindent 6. Datos del sinodal 4 \\ \@supdos \\[1em]	
	\noindent 7. Datos del trabajo escrito \\ \@title \\ \@pags \\ \@degreedate \\[1em]		
	
	\vfill
	\end{alwayssingle}

\fi
}
	

%---------------------------------------------------------------------------------
%                                 Dedication	 						| 
%---------------------------------------------------------------------------------
%The dedication environment makes sure the dedication gets its  own page and is set out in verse format.

\newenvironment{dedication}
{
\begin{alwayssingle}
	\pagestyle{empty}
	\vspace*{1.5cm}
	 {\large{}}
	\vspace{0.5cm}
\begin{flushright}
	\itshape
}
{
\end{flushright}
\end{alwayssingle}
}

%---------------------------------------------------------------------------------
%                                 Acknowledgments 						| 
%---------------------------------------------------------------------------------

	%-------------------------acknowledgements
	% The acknowledgements environment puts a large, bold, centered
	% "Acknowledgements" label at the top of the page. The acknowledgements
	% themselves appear in a quote environment, i.e. tabbed in at both sides, and
	% on its own page.

		\newenvironment{acknowledgements}
		{
		\pagestyle{empty}
		\begin{alwayssingle}
			\vspace*{1.5cm}
			{\chapter*{Agradecimientos}}
			\vspace{0.5cm}
		}
		{\end{alwayssingle}}
	
	%-------------------------acknowledgementslong
	% The acknowledgementslong environment puts a large, bold, centered
	% "Acknowledgements" label at the top of the page. The acknowledgement itself
	% does not appears in a quote environment so you can get more in.

		\newenvironment{acknowledgementslong}
		{
		\pagestyle{empty}
		\begin{alwayssingle}
			\vspace*{1.5cm}
			{\chapter*{Reconocimientos}}
			\vspace{0.5cm}
		}
		{\end{alwayssingle}}

%---------------------------------------------------------------------------------
%                                 Abstract 							| 
%---------------------------------------------------------------------------------

	%-------------------------abstracts
	%The abstract environment puts a large, bold, centered "Abstract" label at
	%the top of the page. The abstract itself appears in a quote environment,
	%i.e. tabbed in at both sides, and on its own page.

		\newenvironment{abstracts}
		{
		\begin{alwayssingle}
			\pagestyle{empty}
			 \begin{center}
		  		\vspace*{1.5cm}
				{\chapter*{Resumen/Abstract}}
			  \end{center}
			  \vspace{0.5cm}
		}
		{\end{alwayssingle}}

	%-------------------------abstractslong
	%The abstractlong environment puts a large, bold, centered "Abstract" label at
	%the top of the page. The abstract itself does not appears in a quote
	%environment so you can get more in.

	\newenvironment{abstractslong} 
	{
	\begin{alwayssingle} 
		\pagestyle{empty}
		\begin{center}
			\vspace*{1.5cm}
			{\chapter*{Resumen}}
		\end{center}
		\vspace{0.5cm}
	}
	{\end{alwayssingle}}

	%------------------------- abstractseparate
	%The abstractseparate environment is for running of a page with the abstract
	%on including title and author etc as required to be handed in separately

	\newenvironment{abstractseparate} {\begin{alwayssingle} \pagestyle{empty}
 	 \vspace*{-1in}
 	\begin{center}
  	  { \Large {\bfseries {\@title}} \par}
  	  {{\large \vspace*{1ex} \@author} \par}
	{\large \vspace*{1ex}
 	   {{UNAM} \par}
  	  {{Facultad de Ciencias} \par}
	\vspace*{1ex}
  	  {{\it \submittedtext} \par}
 	   {\it {\@degree} \par}
	\vspace*{2ex}
  	  {\@degreedate}}
 	 \end{center}}
	{\end{alwayssingle}}


%---------------------------------------------------------------------------------
%                                 Statement of originality if required			| 
%---------------------------------------------------------------------------------

\newenvironment{declaration} {\begin{alwayssingle} \pagestyle{empty}
  \begin{center}
  \vspace*{1.5cm}
  {\chapter*{Declaración de autenticidad}}
  \end{center}
  \vspace{0.5cm}}
{\vspace{1.5cm}\flushright{\@author. \@lugar, \@degreedate}\end{alwayssingle}}


%:-------------------------- page numbers: roman+arabic -----------------------

% ROMANPAGES
%
% The romanpages environment set the page numbering to lowercase roman one
% for the contents and figures lists. It also resets
% page-numbering for the remainder of the dissertation (arabic, starting at 1).

%\newenvironment{romanpages}
%{
%	\setcounter{page}{1}
%	\renewcommand{\thepage}{\roman{page}}
%} % close romanpage env't

{\newpage\renewcommand{\thepage}{\arabic{page}}\setcounter{page}{1}}
