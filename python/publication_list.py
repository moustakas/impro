#!/usr/local/bin/python

"""
Get publication information from ADS and reformat it for LaTeX

Usage:  ./publication_list.py "brammer, g" [--sort-citations --full-header] > mypubs.tex

"""
import sys
import os

import urllib2
import urllib
from xml.dom import minidom

import numpy as np
# import matplotlib.pyplot as plt

def get_citations(author_name='brammer,g', html=False):
    """
    Get XML list of citations for author query, `author_name`, and parse
    it into a format suitable for LaTeX
    """
    ### ADS URL to generate XML, sorted by date
    import locale
    language, output_encoding = locale.getdefaultlocale()
    
    #### The search URL for ADS
    url_template = "http://adsabs.harvard.edu/cgi-bin/nph-abs_connect?db_key=AST&db_key=PRE&qform=AST&arxiv_sel=astro-ph&arxiv_sel=cond-mat&arxiv_sel=cs&arxiv_sel=gr-qc&arxiv_sel=hep-ex&arxiv_sel=hep-lat&arxiv_sel=hep-ph&arxiv_sel=hep-th&arxiv_sel=math&arxiv_sel=math-ph&arxiv_sel=nlin&arxiv_sel=nucl-ex&arxiv_sel=nucl-th&arxiv_sel=physics&arxiv_sel=quant-ph&arxiv_sel=q-bio&sim_query=YES&ned_query=YES&adsobj_query=YES&aut_logic=OR&obj_logic=OR&author=XXXauthorXXX&object=&start_mon=&start_year=&end_mon=&end_year=&ttl_logic=OR&title=&txt_logic=OR&text=&nr_to_return=200&start_nr=1&jou_pick=ALL&article_sel=NO&ref_stems=&data_and=ALL&group_and=ALL&start_entry_day=&start_entry_mon=&start_entry_year=&end_entry_day=&end_entry_mon=&end_entry_year=&min_score=&sort=NDATE&data_type=SHORT_XML&aut_syn=YES&ttl_syn=YES&txt_syn=YES&aut_wt=1.0&obj_wt=1.0&ttl_wt=0.3&txt_wt=3.0&aut_wgt=YES&obj_wgt=YES&ttl_wgt=YES&txt_wgt=YES&ttl_sco=YES&txt_sco=YES&version=1"

    ### Format the search URL for HTML characters
    author_reform = author_name.replace(',','%2C+').replace(' ','+')
    url = url_template.replace('XXXauthorXXX', author_reform)

    if not os.path.exists(author_name.replace(' ','').lower()+'.xml'):

        sys.stderr.write('Retrieve publications for "%s" from http://adsabs.harvard.edu\n' %(author_name))
        
        #### Retrieve the ADS XML query
#       sys.stderr.write(url)
        f = urllib2.urlopen(url)
        xml = f.read()
        f.close()

        #### Seems to be necessary to handle the Unicode characters
        xml = xml.decode(output_encoding)

        #### Translate Unicode characters to LaTeX
        chars = {}
        if html:
            c1, c2 = r'&', r';'
        else:
            c1, c2 = r'$\\', r'$'
            
        chars[u'\u03a9'] = r'%sOmega%s' %(c1, c2)
        chars[u'\u03b1'] = r'%salpha%s' %(c1, c2)
        chars[u'\u03b2'] = r'%sbeta%s' %(c1, c2)
        chars[u'\u03b3'] = r'%sgamma%s' %(c1, c2)
        chars[u'\u03b4'] = r'%sdelta%s' %(c1, c2)
        chars[u'\u03b5'] = r'%sepsilon%s' %(c1, c2)
        chars[u'\u03b6'] = r'%szeta%s' %(c1, c2)
        chars[u'\u03b7'] = r'%seta%s' %(c1, c2)
        chars[u'\u03b8'] = r'%stheta%s' %(c1, c2)
        chars[u'\u03b9'] = r'%siota%s' %(c1, c2)
        chars[u'\u03ba'] = r'%skappa%s' %(c1, c2)
        chars[u'\u03bb'] = r'%slambda%s' %(c1, c2)
        chars[u'\u03bc'] = r'%smu%s' %(c1, c2)
        chars[u'\u03be'] = r'%sxi%s' %(c1, c2)
        chars[u'\u03bf'] = r'%somicron%s' %(c1, c2)
        chars[u'\u03c0'] = r'%spi%s' %(c1, c2)
        chars[u'\u03c1'] = r'%srho%s' %(c1, c2)
        chars[u'\u03c3'] = r'%ssigma%s' %(c1, c2)
        chars[u'\u03a3'] = r'%sSigma%s' %(c1, c2)
        chars[u'\u03c4'] = r'%stau%s' %(c1, c2)
        chars[u'\u03c5'] = r'%snu%s' %(c1, c2)
        chars[u'\u03c6'] = r'%sphi%s' %(c1, c2)
        chars[u'\u03c7'] = r'%schi%s' %(c1, c2)

        chars[u'\xc1'] = r'\'A'
        chars[u'\xe0'] = r'\`a'
        chars[u'\xe1'] = r'\'a'
        chars[u'\xe4'] = r'\"a'
        chars[u'\xe7'] = r'\c{c}'
        chars[u'\xe8'] = r'\`e'
        chars[u'\xe9'] = r'\'e'
        chars[u'\xeb'] = r'\"e'
        chars[u'\xec'] = r'\`i'
        chars[u'\xed'] = r'\'i'
        chars[u'\xef'] = r'\"i'
        chars[u'\xd3'] = r'\'O'
        chars[u'\xd6'] = r'\"O'
        chars[u'\xf1'] = r'\~{n}'
        chars[u'\xf2'] = r'\`o'
        chars[u'\xf3'] = r'\'o'
        chars[u'\xf6'] = r'\"o'
        chars[u'\xf8'] = r'{\o}'
        chars[u'\xf9'] = r'\`u'
        chars[u'\xfa'] = r'\'u'
        chars[u'\xfc'] = r'\"u'
        chars[u'\xc5'] = '\AA\ '
        chars[u'\xcd'] = r'\'I'
        chars[u'\u2272'] = r'x_lt_x'  ### lesssim
        chars[u'\u2273'] = r'x_gt_x'  ### lesssim
        chars[u'\u2a89'] = r'x_lt_x'
        chars[u'\u2a8a'] = r'x_gt_x'
        chars[u'\u2248'] = r'$\approx$'
        chars[u'\xd7'] = r'$\times$'
        chars[u'\u02dc'] = '~'
        chars[u'\u2013'] = '--'
        chars[u'\u201c'] = '``'
        chars[u'\u201d'] = '\'\''
        chars[u'\xdf'] = r'{\ss}'
        chars[u'\u2020'] = r'\dagger'
        chars[u'\xa3'] = r'\pounds\ '
        
        sys.stderr.write('Hey!\n')
        
#       for key in chars.keys():
#           sys.stderr.write('Hey "%s" you\n' %(chars[key]))
#           xml = xml.replace(key, chars[key])
        
        ### Save it to a file.  
        fp = open(author_name.replace(' ','').lower()+'.xml','w')
        fp.write(xml)
        fp.close()
    
    else: 
        ### If the xml file for the input author_name already exists, read it
        fp = open(author_name.replace(' ','').lower()+'.xml','r')
        xml = ''.join(fp.readlines())
        fp.close()

    ### Parse the XML    
    dom = minidom.parseString(xml)
    
    ### Translate bibcode journal names for LaTeX.  
    ### ! Only journals listed here will be added to the output list so 
    ### ! add your own if necessary
    journals = {}
    journals['ApJ'] = 'ApJ'
    journals['ApJS'] = 'ApJS'
    journals['A&A'] = 'A\&A'
    journals['PASP'] = 'PASP'
    journals['MNRAS'] = 'MNRAS'
    journals['Msngr'] = 'The Messenger'
    journals['arXiv1'] = 'arXiv'
    
    #### Some LaTeX characters to replace in the paper titles
    title_chars = {}
    if html:
        title_chars[u'~'] = r'&#126;'
        title_chars[u' > '] = r'&gt;'
        title_chars[u' < '] = r'&lt;'
        title_chars[u'x_lt_x'] = r'&lt;'
        title_chars[u'x_gt_x'] = r'&gt;'
        title_chars[u' <= '] = r'&lt;'
        title_chars[u' >= '] = r'&gt;'
        title_chars[u'z>'] = r'z &gt; '
        title_chars[u'z<'] = r'z &lt; '
        title_chars[u' z '] = r' z '
    else:
        title_chars[u'~'] = r'$\sim$'
        title_chars[u' > '] = r'$>$'
        title_chars[u' < '] = r'$<$'
        title_chars[u'x_lt_x'] = r'$<$'
        title_chars[u'x_gt_x'] = r'$>$'
        title_chars[u' <= '] = r'$<$'
        title_chars[u' >= '] = r'$>$'
        title_chars[u'z>'] = r'$z>$'
        title_chars[u'z<'] = r'$z<$'
        title_chars[u' z '] = r' $z$ '
        
    article_list = []
    cites = []

    #### Loop through the articles in the XML DOM and extract the 
    #### title, bibcode, journal, authors, and number of citations
    articles = dom.getElementsByTagName('record')

#   outfile = open('collaborators.txt','w')
    for article in articles:
        
        ### Title
        try:
            title = article.getElementsByTagName('title')[0].childNodes[0].nodeValue
        except:
            continue
        
        ## Skip if title includes "Erratum"
        if 'Erratum' in title:
            continue
        
        ## Translate some characters in the title to LaTeX format    
        for char in title_chars.keys():
            title = title.replace(char, title_chars[char])
        
        ### Bibcode for URL, also contains Year/Journal/Volume/Page
        try:
            bibcode = article.getElementsByTagName('bibcode')[0].childNodes[0].nodeValue
        except:
            print 'Bibcode failed (%s)' %title
            continue
        
        year = bibcode[0:4]
        journal = bibcode[4:10].replace('.','')
        if journal not in journals.keys():
            print 'Journal not found (%s) in %s' %(journal, title)
            continue
        
        volume = bibcode[10:13].replace('.','')
        page = bibcode[15:18].replace('.','')
        #print bibcode
        
        ### Authors
        try:
            authors = article.getElementsByTagName('author')
        except:
            print 'Authors failed (%s)' %(title)
            pass
        
        ## Make a string list of the author names and make own name bold
        myLast = author_name.split(',')[0].title()
        author_string = ''
        for author in authors:
            auth = author.childNodes[0].nodeValue.split(',')
            last = auth[0].strip()
            first = ','.join(auth[1:]).strip()
            first = first[0]+'.'
            if myLast in last:
                if html:
                    author_string += r'<b>%s %s</b>, ' %(first, last)
                else:
                    author_string += r'\textbf{%s %s}, ' %(first, last)
                    
            else:
                author_string += '%s %s, ' %(first, last)

#               outfile.write('%s %s, ' %(first, last),'\n')

        author_string = author_string[:-2]       
        
        ### Make a citation item for a given article with the LaTeX "authors"
        ### list and "adslink" as defined below in the LaTeX header
        if html:
            test = """
            <p class="pubTitle">
            	"A Strongly-Lensed Massive Ultra-Compact Quiescent Galaxy at $z$ $\sim$ 2.4 in the COSMOS/UltraVISTA Field"
            	<span class="auth"> Muzzin, Labbe, Franx, van Dokkum, Holt, Szomoru, van de Sande, <b>Brammer</b>, Marchesini, Stefanon, Buitrago, Caputi, Dunlop, Fynbo, Le Fevre, McCracken, Milvang-Jensen, 2012, <i>ApJ</i>-submitted (<a href=http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2011arXiv1112.0313W&db_key=PRE&link_type=ABSTRACT>ADS</a>)
            	</span>
            </p>
            """
            
            pub_item = r"""
<p class="pubTitle">
    "%s"
    <span class="auth"> %s, %s, <i>%s</i>, %s, %s (<a href=http://adsabs.harvard.edu/abs/%s>ADS</a>)
    </span>
</p>
            """ %(title, author_string, year, 
                               journals[journal], volume, page, bibcode)
            
        else:
            pub_item = r"""
\item\textit{%s}
\nopagebreak
\begin{authors}
    %s, %s, \textit{%s}, %s, %s
    \adslink{http://adsabs.harvard.edu/abs/%s}
\end{authors}""" %(title, author_string, year, 
                   journals[journal], volume, page, bibcode)
        
        ### Get number of citations
        try:
            cites.append( article.getElementsByTagName('citations')[0].childNodes[0].nodeValue)
        except:
            cites.append(0)
        
        ### Add the article to the list
        article_list.append(pub_item)

#   outfile.close()
    
    ### cites -> integer
    cites = np.cast[int](cites)

    ### Done!
    return article_list, cites
    
def run_it(argv):
    """
    Retrieve the article list and print the LaTeX-formatted list to 
    the terminal
    """
    ### Usage test
    if len(argv) == 1:
        print '\nUsage:  ./publication_list.py "brammer, g" [--sort-citations --full-header --html] > mypubs.tex\n'
        return False
    
    #### sort the output by citations and add the number of citations to the 
    #### title string    
    if (len(argv) > 2) & ('--sort-citations' in argv):
        sort_cite = True
    else:
        sort_cite = False
    
    #### HTML formatting rather than Latex
    if (len(argv) > 2) & ('--html' in argv):
        format_html = True
    else:
        format_html = False
    
    #### Make a full stand-alone LaTeX file that can be compiled itself.  
    #### Otherwise, just include the publication list itself.
    if (len(argv) > 2) & ('--full-header' in argv):
        full_head = True
    else:
        full_head = False
    
    ### Get the list here
    try:
        article_list, ncites = get_citations(author_name=argv[1], html=format_html)
    except:
        sys.stderr.write("Encountered an error retrieving publication list for \"%s\"\n" %(argv[1]))
        os.remove(argv[1].replace(' ','')+'.xml')
        return False
    
    #### Simple HTML format, note latex special characters not translated yet
    if format_html:
        if sort_cite:
            so = np.argsort(ncites)[::-1]
        else:
            so = np.arange(len(article_list))

        ### Print the article information
        for i in so:
            if sort_cite:
                line = article_list[i].replace('pubTitle">', 'pubTitle"><b>%d</b>' %(ncites[i]))
            else:
                line = article_list[i]
            #
            print line
        
        return True
        
    ### Full stand-alone header    
    full_header = r"""

\documentclass[12pt]{report}

\usepackage{graphicx}

\usepackage{color,hyperref}
\definecolor{YaleBlue}{rgb}{0.0,0.22,0.444}

\hypersetup{colorlinks,breaklinks,
            linkcolor=YaleBlue,urlcolor=YaleBlue,
            anchorcolor=YaleBlue,citecolor=YaleBlue}

\usepackage{enumitem}

\usepackage{multicol}

%%%% Use sans-serif fonts
\renewcommand{\familydefault}{\sfdefault}

\addtolength{\textwidth}{1in}
\addtolength{\hoffset}{-0.5in}

% 
\addtolength{\textheight}{4cm}
\addtolength{\voffset}{-2cm}

\newlength{\secoffset}\setlength{\secoffset}{8pt}%

\reversemarginpar

\setlength{\parindent}{0in}

\thispagestyle{empty}

\renewcommand{\section}[1]%
{\phantomsection\addcontentsline{toc}{subsection}{#1}%
 \vspace{0.5cm}
% {\large { \scshape #1} }
{ \scshape #1}

 \vspace{0.2cm}
}

\begin{document} 
"""
    
    ### Just the information needed for the "authors" and "my_enumerate" 
    ### functions
    header = r"""\newcommand{\firstauthor}
		{\pagebreak[2]
 		\hspace{0in}%
 		\marginpar{\raggedleft \ding{192}}}


\newenvironment{my_enumerate}{
	\begin{enumerate}[leftmargin=0pt,topsep=0pt,noitemsep]
	}{\end{enumerate}}
\newcommand\litem[0]{\item[\alph{enumi}]}

%\usepackage[T1]{fontenc}
%\usepackage[scaled]{helvet}
%\renewcommand*\familydefault{\sfdefault}

% To add some paragraph space between lines.
% This also tells LaTeX to preferably break a page on one of these gaps
% if there is a needed pagebreak nearby.
\newcommand{\blankline}{\quad\pagebreak[2]}

% Author list for Publications
\newenvironment{authors}%
{\nopagebreak %\hspace{-10pt}
\begin{itemize}[topsep=0pt]{}%
         \item[]{\setlength{\leftmargin}{1cm}}%
         %
}
{\end{itemize}\smallskip}

\newcommand{\adslink}[1]{\href{#1}{\scriptsize{\textsf{[ADS]}}}}

\section{Publications} 

\begin{my_enumerate}"""
    
    footer = r"""
\end{my_enumerate}"""

    ### Sorting by citations
    if sort_cite:
        so = np.argsort(ncites)[::-1]
    else:
        so = np.arange(len(article_list))
    
    ### Print the full header
    if full_head:
        print full_header

    print header
    
    ### Print the article information
    for i in so:
        if sort_cite:
            line = article_list[i].replace(r'\item\textit{',r'\item{\bf {\tiny [%d]\ }}\textit{' %(ncites[i]))
        else:
            line = article_list[i]

        print line
    
    ### LaTex end things
    print footer   
    if full_head:
        print r"""
\end{document}
"""

if __name__ == "__main__":
    run_it(sys.argv)
    
