tugHall version 3.0
====================
[![License](https://img.shields.io/badge/License-GPLv3-orange.svg)](https://github.com/tugHall/tugHall.3/blob/main/LICENSE.md)

**tugHall** _(**tu**mor **g**ene-**Hall**mark)_ is a cancer-cell evolution model simulator, wherein gene mutations are linked to the hallmarks of cancer, 
which influence tumor cell behaviors. 

This is an _**R**_-based package to simulate the cancer cell evolution in the framework of the model proposed by _**Prof. Mamoru Kato**_,
_Head of Bioinformatics Division, Research Institute, National Cancer Center, Tokyo, JAPAN_.

Authors and contributor list:
---
_**Iurii Nagornov**_ (Maintainer, Author)

_**Munmee Dutta**_ (Contributor)

_**Jo Nishino**_ (Contributor)

_**Eisaku Furukawa**_ (Contributor)

_**Mamoru Kato**_ (Author)

_Division of Bioinformatics, Research Institute, National Cancer Center Japan, Tokyo, Japan_

All questions and requests can be sent to inagonov@ncc.go.jp or nagornov.yuri@gmail.com 

Short description
---
The wide availability of recent cancer genomic data requires a coherent model that can sort out the relevant findings to systematically explain the clonal evolution and resultant intra-tumor heterogeneity (ITH). Here, we present a new mathematical model designed to computationally simulate the evolution of cancer cells. The model connects well-known cancer hallmarks with the specific mutational states of tumor-related genes. The cell behavior phenotypes are stochastically determined and the hallmarks interfere probabilistically with the phenotypic probabilities. In turn, the hallmark variables depend on the mutational states of tumor-related genes. Thus, it is expected our software can be used to deepen our understanding of cancer-cell evolution and generation of ITH.

Installation
---
Please, kindly install tugHall.3 package from source file.
Firstly, download file 'tugHall.3_3.0.tar.gz' from repository to the working folder, and then:

    install.packages("./tugHall.3_3.0.tar.gz", repos = NULL, type = "source")

How to use **tugHall.3** and how to analyze data, kindly see vignettes in the 
package documentation or see code in the pipeline examples. To get pipelines please, run the function:

    copy_pipelines( dir = './')

## Cite package tugHall.3

For publication, please, be kind to use next references related to tugHall software:

- Model description and first version of tugHall

Iurii S Nagornov,  Mamoru Kato. tugHall: a simulator of cancer-cell evolution based on the hallmarks of cancer and tumor-related genes. 
[Bioinformatics, V.36, N11, June 2020, pp. 3597–3599](https://doi.org/10.1093/bioinformatics/btaa182)

- Clone-based version of tugHall 2.0 and 2.1 with acceleration of calculation speed

Nagornov, I., Nishino, J., Kato, M. (2020). tugHall: A Tool to Reproduce Darwinian Evolution of Cancer Cells for Simulation-Based Personalized Medicine. In: Mathematical and Computational Oncology. ISMCO 2020. 
[Lecture Notes in Computer Science, vol 12508. Springer](https://doi.org/10.1007/978-3-030-64511-3_7)

- for reference of tugHall.3 package, please, cite both papers and corresponding github repository: 

https://github.com/tugHall/tugHall.3 


