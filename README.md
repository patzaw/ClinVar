-   [General description](#general-description)
-   [DESCRIPTION](#description)
-   [Data model](#data-model)
-   [Scripts](#scripts)
-   [DODO](#dodo)

General description
===================

This repository contains scripts used to parse the
*ClinVarFullRelease.xml* file from
[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) and organize the data
according to a customized data model.

> This is **not** an official repository for ClinVar, please refer to
> the website of [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/).

DESCRIPTION
===========

The **DESCRIPTION.json** file provides a formal description of the
project including its version, its maintainer and the original sources
of information.

Data model
==========

The data model has been formalized using the
[ReDaMoR](https://github.com/patzaw/ReDaMoR) R package (available on
[CRAN](https://cran.r-project.org/package=ReDaMoR)). It is available in
the **model** folder.

The **Collections** subfolder contains json files which formally
describe tables gathering information about key concepts, such as genes
or diseases, that can be used create cross-reference with other
resources. **Collections** mechanisms are developed in the frame of the
[TKCat](https://github.com/patzaw/TKCat) project.

Scripts
=======

Parsing scripts are located in the **scripts** folder.

DODO
====

The created flat files are used to feed the database supporting DODO
([Dictionary of Disease Ontologies](https://github.com/Elysheba/DODO)).
