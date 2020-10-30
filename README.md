# Graphical-Set-based model
My thessis on creating a new information retrieval model representing textual documents as graphs, considering node importance, the retrieval proccess is handled by the set based model. For more information:  https://rdcu.be/b8Y1S

This is the code for the paper "**A Graph-Based Extension for the Set-Based Model Implementing Algorithms Based on Important Nodes**" presented on _Artificial Intelligence Applications and Innovations, AIAI 2020 IFIP WG 12.5 International Workshops_.

## Implementation - Description

At this point, we dive into the expirimental proccess. For detailed algorithmic explanation review the above mentioned paper. The testing proccess consist of 4 indexing parts and one query, functions declared with a simple 1-2 menu option. At first, through a menu function we declare which models the indexing proccess would like to implement. _Important note: For each indexing proccess the importance of keywords (h) and percentage of prunning (P) parameters should be the same._ 

After the 4 inverted indexes are created through the starting menu function we declare the start of the query proccess, which loads the 4 inverted indexes and calulate the respective weights, as well as the simple set-based model weights, creating thus creating 5 comparison models, having set-based model as a baseline. The results are stored in an .xls file _(Important node: path files and file names should be fixed.)_ The results are focused on recall and precision of each query out of 100. A simple graphical approach was used too.        

## Dataset

As expirimental dataset the Cystic Fibrosis collection was used (REF: http://www2.dcc.ufmg.br/livros/irbook/cfc.html ). The CF collection consists of 1239 documents-articles, including author, references, citations, abstract and more information. Also it includes 100 queries and the most relevant documents for each query. The ranking proccess handled by experts using "a grade" ranking proccess.

## How to use it

At first we need to extract the main text from our texual data. The main text is expressed by the abstract of each document, thus the "tesxtparser.py" is used. Using specific acronyms which the CF collection supports it extracts each part of each document and stores the needed data in a specific path. _IMPORTANT note: paths need to be altered depending on your machine_.

After, the "preproccess.py" handles the data preproccess functions. In this step, the tokenization and stopword filltering  proccesses are taking place.
The queries is extracted with "queryparser.py" and pasted as a list in the "main.py"-main function. The main function was used through the mainautomation script which inputs as arguments the testing parameters. 

There are **two** key points that require attention. Firstly, the inverted indexes should get renamed after each expiriment because the new inverted index will not overwrite the old one but it will get appended in the old file. Last but not least, on each collection change, such as different preproccess the docinfo.dat should get manually cleared or deleted. 

## Citation

Please cite as:

_Simple citation:_

        Kalogeropoulos NR., Doukas I., Makris C., Kanavos A. (2020) A Graph-Based Extension for the Set-Based Model Implementing
        Algorithms Based on Important Nodes. In: Maglogiannis I., Iliadis L., Pimenidis E. (eds) Artificial Intelligence 
        Applications and Innovations. AIAI 2020 IFIP WG 12.5 International Workshops. AIAI 2020. IFIP Advances in Information
        and Communication Technology, vol 585. Springer, Cham. https://doi.org/10.1007/978-3-030-49190-1_13

_Bibtex citation:_

    @InProceedings{10.1007/978-3-030-49190-1_13,
    author="Kalogeropoulos, Nikitas-Rigas
            and Doukas, Ioannis
            and Makris, Christos
            and Kanavos, Andreas",
            editor="Maglogiannis, Ilias
            and Iliadis, Lazaros
            and Pimenidis, Elias",
    title="A Graph-Based Extension for the Set-Based Model Implementing Algorithms Based on Important Nodes",
    booktitle="Artificial Intelligence Applications and Innovations. AIAI 2020 IFIP WG 12.5 International Workshops",
    year="2020",
    publisher="Springer International Publishing",
    address="Cham",
    pages="143--154",
    isbn="978-3-030-49190-1"
     }



