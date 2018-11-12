Cancer predisposition genes
---------------------------

In order to establish a list of cancer predisposition genes, we considered both established, manually curated sources, and also empirical data (pathogenic gene variants) associated with hereditary cancer phenotypes.

Initially, we used three manually curated sources to compile a list of 209 protein-coding genes with associations to cancer predisposition/cancer syndromes:

-   *TCGA\_PANCAN\_18* - TCGCA Pancancer germline study - [Huang et al, Cell, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052)
-   *CGC\_86* - Curated list - [Cancer Gene Census (COSMIC)](https://cancer.sanger.ac.uk/census) - version 86
-   *NCGC* - Expert-curated list from [Norwegian Cancer Genomics Consortium](http://cancergenomics.no)

Data with respect to mechanisms of inheritance (<i>MoI</i> - autosomal recessive (AR) vs. autosomal dominant (AD) etc.) and whether mechanisms of disease are associated with loss-of-function (<i>LoF</i>) or gain-of-function (<i>GoF</i>) were primarily retrieved from [Maxwell et al., Am J Hum Genet, 2016](https://www.ncbi.nlm.nih.gov/pubmed/27153395).

<!-- Next, we considered the tree of [Medical Subject Headings (MeSH)](https://www.nlm.nih.gov/mesh/intro_trees.html) in [MedGen](https://www.ncbi.nlm.nih.gov/medgen/), and curated a list of hereditary cancer phenotypes, using *[Inherited cancer-predisposing syndrome](https://www.ncbi.nlm.nih.gov/medgen/798871)* and *[Hereditary Cancer](https://www.ncbi.nlm.nih.gov/medgen/232504)* as the main starting points. For this curated list of inherited cancer phenotypes, we queried ClinVar for protein-coding genes with class4/class5 variants. This resulted in a list of 166 protein-coding genes. Finally, limiting this to the ones that were also listed in the non-redundant set from three manually curated sources, our final list contained 144 protein-coding genes. -->
We want to make it explicit that this list of 209 genes is by no means regarded as an international consensus, but should rather be subject to continuous update by the international community that carry expertise on genetic risk factors for cancer.

<table style="width:56%;">
<colgroup>
<col width="6%" />
<col width="6%" />
<col width="6%" />
<col width="6%" />
<col width="6%" />
<col width="6%" />
<col width="6%" />
<col width="6%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Symbol</th>
<th align="left">Entrez ID</th>
<th align="left">MOI</th>
<th align="left">LoF</th>
<th align="left">Gene Name</th>
<th align="left">Source</th>
<th align="left">Phenotype_Syndrome_CUI</th>
<th align="left">Phenotype_Syndrome_Term</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">ABCB11</td>
<td align="left">8647</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">ATP binding cassette subfamily B member 11</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C2931132</td>
<td align="left">Crigler-Najjar syndrome, type II</td>
</tr>
<tr class="even">
<td align="left">ABRAXAS1</td>
<td align="left">84142</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">abraxas 1, BRCA1 A complex subunit</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">ACD</td>
<td align="left">65057</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">ACD, shelterin complex subunit and telomerase recruitment factor</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">AIP</td>
<td align="left">9049</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">aryl hydrocarbon receptor interacting protein</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">AKT1</td>
<td align="left">207</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">AKT serine/threonine kinase 1</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">ALK</td>
<td align="left">238</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">ALK receptor tyrosine kinase</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0027819</td>
<td align="left">Neuroblastoma</td>
</tr>
<tr class="odd">
<td align="left">APC</td>
<td align="left">324</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">APC, WNT signaling pathway regulator</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C2713442</td>
<td align="left">Familial adenomatous polyposis 1</td>
</tr>
<tr class="even">
<td align="left">APOBEC3B</td>
<td align="left">9582</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">apolipoprotein B mRNA editing enzyme catalytic subunit 3B</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">AR</td>
<td align="left">367</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">androgen receptor</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">ATM</td>
<td align="left">472</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">ATM serine/threonine kinase</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0004135</td>
<td align="left">Ataxia-telangiectasia syndrome</td>
</tr>
<tr class="odd">
<td align="left">ATR</td>
<td align="left">545</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">ATR serine/threonine kinase</td>
<td align="left">TCGA_PANCAN_2018,CGC_86</td>
<td align="left">C0265202;C3281203</td>
<td align="left">Seckel syndrome; Cutaneous telangiectasia and cancer syndrome, familial</td>
</tr>
<tr class="even">
<td align="left">AXIN1</td>
<td align="left">8312</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">axin 1</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">AXIN2</td>
<td align="left">8313</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">axin 2</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1837750</td>
<td align="left">Oligodontia-colorectal cancer syndrome</td>
</tr>
<tr class="even">
<td align="left">BAP1</td>
<td align="left">8314</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">BRCA1 associated protein 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C3280492;CN235077</td>
<td align="left">Tumor susceptibility linked to germline BAP1 mutations; BAP1 Cancer Syndrome</td>
</tr>
<tr class="odd">
<td align="left">BARD1</td>
<td align="left">580</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">BRCA1 associated RING domain 1</td>
<td align="left">CGC_86,NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">BLM</td>
<td align="left">641</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">Bloom syndrome RecQ like helicase</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0005859</td>
<td align="left">Bloom syndrome</td>
</tr>
<tr class="odd">
<td align="left">BMPR1A</td>
<td align="left">657</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">bone morphogenetic protein receptor type 1A</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0345893</td>
<td align="left">Juvenile polyposis syndrome</td>
</tr>
<tr class="even">
<td align="left">BRAF</td>
<td align="left">673</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">B-Raf proto-oncogene, serine/threonine kinase</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C3150970;C3150971</td>
<td align="left">Noonan syndrome 7; LEOPARD syndrome 3</td>
</tr>
<tr class="odd">
<td align="left">BRCA1</td>
<td align="left">672</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">BRCA1, DNA repair associated</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0677776</td>
<td align="left">Hereditary breast and ovarian cancer syndrome</td>
</tr>
<tr class="even">
<td align="left">BRCA2</td>
<td align="left">675</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">BRCA2, DNA repair associated</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0677776;C1838457</td>
<td align="left">Hereditary breast and ovarian cancer syndrome; Fanconi anemia, complementation group D1</td>
</tr>
<tr class="odd">
<td align="left">BRIP1</td>
<td align="left">83990</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">BRCA1 interacting protein C-terminal helicase 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1836860</td>
<td align="left">Fanconi anemia, complementation group J</td>
</tr>
<tr class="even">
<td align="left">BUB1B</td>
<td align="left">701</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">BUB1 mitotic checkpoint serine/threonine kinase B</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1850343</td>
<td align="left">Mosaic variegated aneuploidy syndrome</td>
</tr>
<tr class="odd">
<td align="left">CASR</td>
<td align="left">846</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">calcium sensing receptor</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">CBL</td>
<td align="left">867</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">Cbl proto-oncogene</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C0028326</td>
<td align="left">Noonan syndrome</td>
</tr>
<tr class="odd">
<td align="left">CDC73</td>
<td align="left">79577</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">cell division cycle 73</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1704981</td>
<td align="left">Hyperparathyroidism 2</td>
</tr>
<tr class="even">
<td align="left">CDH1</td>
<td align="left">999</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">cadherin 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1708349</td>
<td align="left">Hereditary diffuse gastric cancer</td>
</tr>
<tr class="odd">
<td align="left">CDH10</td>
<td align="left">1008</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">cadherin 10</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">CDK4</td>
<td align="left">1019</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">cyclin dependent kinase 4</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1512419;C1836892</td>
<td align="left">Hereditary cutaneous melanoma; Cutaneous malignant melanoma 3</td>
</tr>
<tr class="odd">
<td align="left">CDKN1B</td>
<td align="left">1027</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">cyclin dependent kinase inhibitor 1B</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1970712</td>
<td align="left">Multiple endocrine neoplasia, type 4</td>
</tr>
<tr class="even">
<td align="left">CDKN1C</td>
<td align="left">1028</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">cyclin dependent kinase inhibitor 1C</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C0004903</td>
<td align="left">Beckwith-Wiedemann syndrome</td>
</tr>
<tr class="odd">
<td align="left">CDKN2A</td>
<td align="left">1029</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">cyclin dependent kinase inhibitor 2A</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1512419</td>
<td align="left">Hereditary cutaneous melanoma</td>
</tr>
<tr class="even">
<td align="left">CEBPA</td>
<td align="left">1050</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">CCAAT enhancer binding protein alpha</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C0023467</td>
<td align="left">Acute myeloid leukemia</td>
</tr>
<tr class="odd">
<td align="left">CEP57</td>
<td align="left">9702</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">centrosomal protein 57</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">CHEK2</td>
<td align="left">11200</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">checkpoint kinase 2</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0346153</td>
<td align="left">Familial cancer of breast</td>
</tr>
<tr class="odd">
<td align="left">COL7A1</td>
<td align="left">1294</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">collagen type VII alpha 1 chain</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C0014527</td>
<td align="left">Epidermolysis bullosa</td>
</tr>
<tr class="even">
<td align="left">CTNNA1</td>
<td align="left">1495</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">catenin alpha 1</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">CTNNB1</td>
<td align="left">1499</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">catenin beta 1</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">CTR9</td>
<td align="left">9646</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">CTR9 homolog, Paf1/RNA polymerase II complex component</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">CTRC</td>
<td align="left">11330</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">chymotrypsin C</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">CXCR4</td>
<td align="left">7852</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">C-X-C motif chemokine receptor 4</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">CYLD</td>
<td align="left">1540</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">CYLD lysine 63 deubiquitinase</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1857941</td>
<td align="left">Spiegler-Brooke syndrome</td>
</tr>
<tr class="even">
<td align="left">DDB2</td>
<td align="left">1643</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">damage specific DNA binding protein 2</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1848411</td>
<td align="left">Xeroderma pigmentosum, group E</td>
</tr>
<tr class="odd">
<td align="left">DICER1</td>
<td align="left">23405</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">dicer 1, ribonuclease III</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1266144;C3839822</td>
<td align="left">Pleuropulmonary blastoma; DICER1 syndrome</td>
</tr>
<tr class="even">
<td align="left">DIRAS3</td>
<td align="left">9077</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">DIRAS family GTPase 3</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">DIS3L2</td>
<td align="left">129563</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">DIS3 like 3'-5' exoribonuclease 2</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C0796113</td>
<td align="left">Renal hamartomas nephroblastomatosis and fetal gigantism</td>
</tr>
<tr class="even">
<td align="left">DKC1</td>
<td align="left">1736</td>
<td align="left">XLR</td>
<td align="left">LoF</td>
<td align="left">dyskerin pseudouridine synthase 1</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C1148551</td>
<td align="left">Dyskeratosis congenita X-linked</td>
</tr>
<tr class="odd">
<td align="left">DOCK8</td>
<td align="left">81704</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">dedicator of cytokinesis 8</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C1968689</td>
<td align="left">Hyperimmunoglobulin E recurrent infection syndrome, autosomal recessive</td>
</tr>
<tr class="even">
<td align="left">DROSHA</td>
<td align="left">29102</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">drosha ribonuclease III</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">DTX3L</td>
<td align="left">151636</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">deltex E3 ubiquitin ligase 3L</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">EGFR</td>
<td align="left">1956</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">epidermal growth factor receptor</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">ELANE</td>
<td align="left">1991</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">elastase, neutrophil expressed</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C1859966</td>
<td align="left">Severe congenital neutropenia autosomal dominant</td>
</tr>
<tr class="even">
<td align="left">ENG</td>
<td align="left">2022</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">endoglin</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">EPCAM</td>
<td align="left">4072</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">epithelial cell adhesion molecule</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C2750471</td>
<td align="left">Hereditary nonpolyposis colorectal cancer type 8</td>
</tr>
<tr class="even">
<td align="left">ERBB4</td>
<td align="left">2066</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">erb-b2 receptor tyrosine kinase 4</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">ERCC1</td>
<td align="left">2067</td>
<td align="left">AR</td>
<td align="left">NA</td>
<td align="left">ERCC excision repair 1, endonuclease non-catalytic subunit</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C1853100</td>
<td align="left">Cerebrooculofacioskeletal syndrome 4</td>
</tr>
<tr class="even">
<td align="left">ERCC2</td>
<td align="left">2068</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">ERCC excision repair 2, TFIIH core complex helicase subunit</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0268138</td>
<td align="left">Xeroderma pigmentosum, group D</td>
</tr>
<tr class="odd">
<td align="left">ERCC3</td>
<td align="left">2071</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">ERCC excision repair 3, TFIIH core complex helicase subunit</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0268136</td>
<td align="left">Xeroderma Pigmentosum, Complementation Group B</td>
</tr>
<tr class="even">
<td align="left">ERCC4</td>
<td align="left">2072</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">ERCC excision repair 4, endonuclease catalytic subunit</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0268140;C3808988</td>
<td align="left">Xeroderma pigmentosum, group F; Fanconi anemia, complementation group Q</td>
</tr>
<tr class="odd">
<td align="left">ERCC5</td>
<td align="left">2073</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">ERCC excision repair 5, endonuclease</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0268141</td>
<td align="left">Xeroderma pigmentosum, group G</td>
</tr>
<tr class="even">
<td align="left">ETV6</td>
<td align="left">2120</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">ETS variant 6</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C4015537</td>
<td align="left">Thrombocytopenia 5</td>
</tr>
<tr class="odd">
<td align="left">EXT1</td>
<td align="left">2131</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">exostosin glycosyltransferase 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0008479</td>
<td align="left">Chondrosarcoma</td>
</tr>
<tr class="even">
<td align="left">EXT2</td>
<td align="left">2132</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">exostosin glycosyltransferase 2</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1851413</td>
<td align="left">Multiple exostoses type 2</td>
</tr>
<tr class="odd">
<td align="left">EZH2</td>
<td align="left">2146</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">enhancer of zeste 2 polycomb repressive complex 2 subunit</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">FAH</td>
<td align="left">2184</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">fumarylacetoacetate hydrolase</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">FANCA</td>
<td align="left">2175</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">Fanconi anemia complementation group A</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C3469521</td>
<td align="left">Fanconi anemia, complementation group A</td>
</tr>
<tr class="even">
<td align="left">FANCB</td>
<td align="left">2187</td>
<td align="left">XLR</td>
<td align="left">LoF</td>
<td align="left">Fanconi anemia complementation group B</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">FANCC</td>
<td align="left">2176</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">Fanconi anemia complementation group C</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C3468041</td>
<td align="left">Fanconi anemia, complementation group C</td>
</tr>
<tr class="even">
<td align="left">FANCD2</td>
<td align="left">2177</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">Fanconi anemia complementation group D2</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C3160738</td>
<td align="left">Fanconi anemia, complementation group D2</td>
</tr>
<tr class="odd">
<td align="left">FANCE</td>
<td align="left">2178</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">Fanconi anemia complementation group E</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C3160739</td>
<td align="left">Fanconi anemia, complementation group E</td>
</tr>
<tr class="even">
<td align="left">FANCF</td>
<td align="left">2188</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">Fanconi anemia complementation group F</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C3469526</td>
<td align="left">Fanconi anemia, complementation group F</td>
</tr>
<tr class="odd">
<td align="left">FANCG</td>
<td align="left">2189</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">Fanconi anemia complementation group G</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C3469527</td>
<td align="left">Fanconi anemia, complementation group G</td>
</tr>
<tr class="even">
<td align="left">FANCI</td>
<td align="left">55215</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">Fanconi anemia complementation group I</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C1836861</td>
<td align="left">Fanconi anemia, complementation group I</td>
</tr>
<tr class="odd">
<td align="left">FANCL</td>
<td align="left">55120</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">Fanconi anemia complementation group L</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C3469528</td>
<td align="left">Fanconi anemia, complementation group L</td>
</tr>
<tr class="even">
<td align="left">FANCM</td>
<td align="left">57697</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">Fanconi anemia complementation group M</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">FAS</td>
<td align="left">355</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">Fas cell surface death receptor</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C1328840</td>
<td align="left">Autoimmune lymphoproliferative syndrome</td>
</tr>
<tr class="even">
<td align="left">FAT1</td>
<td align="left">2195</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">FAT atypical cadherin 1</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">FEN1</td>
<td align="left">2237</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">flap structure-specific endonuclease 1</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">FH</td>
<td align="left">2271</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">fumarate hydratase</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1708350</td>
<td align="left">Multiple cutaneous leiomyomas</td>
</tr>
<tr class="odd">
<td align="left">FLCN</td>
<td align="left">201163</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">folliculin</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">CN221571</td>
<td align="left">Birt-Hogg-Dub syndrome</td>
</tr>
<tr class="even">
<td align="left">GALNT12</td>
<td align="left">79695</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">polypeptide N-acetylgalactosaminyltransferase 12</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">GATA2</td>
<td align="left">2624</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">GATA binding protein 2</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C3280030</td>
<td align="left">Dendritic cell, monocyte, B lymphocyte, and natural killer lymphocyte deficiency</td>
</tr>
<tr class="even">
<td align="left">GBA</td>
<td align="left">2629</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">glucosylceramidase beta</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">GJB2</td>
<td align="left">2706</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">gap junction protein beta 2</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">GPC3</td>
<td align="left">2719</td>
<td align="left">XLR</td>
<td align="left">LoF</td>
<td align="left">glypican 3</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0796154</td>
<td align="left">Simpson-Golabi-Behmel syndrome</td>
</tr>
<tr class="odd">
<td align="left">GREM1</td>
<td align="left">26585</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">gremlin 1, DAN family BMP antagonist</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">HABP2</td>
<td align="left">3026</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">hyaluronan binding protein 2</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">HFE</td>
<td align="left">3077</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">homeostatic iron regulator</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C0392514</td>
<td align="left">Hereditary hemochromatosis</td>
</tr>
<tr class="even">
<td align="left">HMBS</td>
<td align="left">3145</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">hydroxymethylbilane synthase</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">HNF1A</td>
<td align="left">6927</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">HNF1 homeobox A</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1840646</td>
<td align="left">Hepatic adenomas, familial</td>
</tr>
<tr class="even">
<td align="left">HNF1B</td>
<td align="left">6928</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">HNF1 homeobox B</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">HOXB13</td>
<td align="left">10481</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">homeobox B13</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">HRAS</td>
<td align="left">3265</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">HRas proto-oncogene, GTPase</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0587248</td>
<td align="left">Costello syndrome</td>
</tr>
<tr class="odd">
<td align="left">ITK</td>
<td align="left">3702</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">IL2 inducible T cell kinase</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C3552634</td>
<td align="left">Lymphoproliferative syndrome 1</td>
</tr>
<tr class="even">
<td align="left">JMJD1C</td>
<td align="left">221037</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">jumonji domain containing 1C</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">KDR</td>
<td align="left">3791</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">kinase insert domain receptor</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">KIF1B</td>
<td align="left">23095</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">kinesin family member 1B</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">KIT</td>
<td align="left">3815</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">KIT proto-oncogene receptor tyrosine kinase</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0238198</td>
<td align="left">Gastrointestinal stroma tumor</td>
</tr>
<tr class="even">
<td align="left">KRAS</td>
<td align="left">3845</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">KRAS proto-oncogene, GTPase</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C1860991;C2674723</td>
<td align="left">Noonan syndrome 3; RAS-associated autoimmune leukoproliferative disorder</td>
</tr>
<tr class="odd">
<td align="left">LMO1</td>
<td align="left">4004</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">LIM domain only 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86</td>
<td align="left">C4225207</td>
<td align="left">Neuroblastoma, susceptibility to, 7</td>
</tr>
<tr class="even">
<td align="left">LZTR1</td>
<td align="left">8216</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">leucine zipper like transcription regulator 1</td>
<td align="left">CGC_86,NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">MAP2K1</td>
<td align="left">5604</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">mitogen-activated protein kinase kinase 1</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">MAP2K2</td>
<td align="left">5605</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">mitogen-activated protein kinase kinase 2</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">MAX</td>
<td align="left">4149</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">MYC associated factor X</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1708353</td>
<td align="left">Hereditary Paraganglioma-Pheochromocytoma Syndromes</td>
</tr>
<tr class="even">
<td align="left">MEN1</td>
<td align="left">4221</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">menin 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0025267</td>
<td align="left">Multiple endocrine neoplasia, type 1</td>
</tr>
<tr class="odd">
<td align="left">MET</td>
<td align="left">4233</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">MET proto-oncogene, receptor tyrosine kinase</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C0007134</td>
<td align="left">Renal cell carcinoma, papillary, 1</td>
</tr>
<tr class="even">
<td align="left">MITF</td>
<td align="left">4286</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">melanogenesis associated transcription factor</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C1512419</td>
<td align="left">Hereditary cutaneous melanoma</td>
</tr>
<tr class="odd">
<td align="left">MLH1</td>
<td align="left">4292</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">mutL homolog 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0009405;C0265325;C1333990</td>
<td align="left">Hereditary nonpolyposis colon cancer; Turcot syndrome; Lynch syndrome</td>
</tr>
<tr class="even">
<td align="left">MLH3</td>
<td align="left">27030</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">mutL homolog 3</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">MPL</td>
<td align="left">4352</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">MPL proto-oncogene, thrombopoietin receptor</td>
<td align="left">TCGA_PANCAN_2018,CGC_86</td>
<td align="left">C4273671</td>
<td align="left">Inherited predisposition to essential thrombocythemia</td>
</tr>
<tr class="even">
<td align="left">MRE11</td>
<td align="left">4361</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">MRE11 homolog, double strand break repair nuclease</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">MSH2</td>
<td align="left">4436</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">mutS homolog 2</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0009405;C0265325;C1333990</td>
<td align="left">Hereditary nonpolyposis colon cancer; Turcot syndrome; Lynch syndrome</td>
</tr>
<tr class="even">
<td align="left">MSH3</td>
<td align="left">4437</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">mutS homolog 3</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">MSH6</td>
<td align="left">2956</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">mutS homolog 6</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0009405;C0265325;C1333990</td>
<td align="left">Hereditary nonpolyposis colon cancer; Turcot syndrome; Lynch syndrome</td>
</tr>
<tr class="even">
<td align="left">MTAP</td>
<td align="left">4507</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">methylthioadenosine phosphorylase</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C1862177</td>
<td align="left">Diaphyseal medullary stenosis with malignant fibrous histiocytoma</td>
</tr>
<tr class="odd">
<td align="left">MUTYH</td>
<td align="left">4595</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">mutY DNA glycosylase</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1837991</td>
<td align="left">MYH-associated polyposis</td>
</tr>
<tr class="even">
<td align="left">NBN</td>
<td align="left">4683</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">nibrin</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0398791</td>
<td align="left">Microcephaly, normal intelligence and immunodeficiency</td>
</tr>
<tr class="odd">
<td align="left">NF1</td>
<td align="left">4763</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">neurofibromin 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0027831</td>
<td align="left">Neurofibromatosis, type 1</td>
</tr>
<tr class="even">
<td align="left">NF2</td>
<td align="left">4771</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">neurofibromin 2</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0027832</td>
<td align="left">Neurofibromatosis, type 2</td>
</tr>
<tr class="odd">
<td align="left">NHP2</td>
<td align="left">55651</td>
<td align="left">AR</td>
<td align="left">NA</td>
<td align="left">NHP2 ribonucleoprotein</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C3151441</td>
<td align="left">Dyskeratosis congenita, autosomal recessive 2</td>
</tr>
<tr class="even">
<td align="left">NOP10</td>
<td align="left">55505</td>
<td align="left">AR</td>
<td align="left">NA</td>
<td align="left">NOP10 ribonucleoprotein</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C1857144</td>
<td align="left">Dyskeratosis congenita autosomal recessive 1</td>
</tr>
<tr class="odd">
<td align="left">NRAS</td>
<td align="left">4893</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">NRAS proto-oncogene, GTPase</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C0265318;C0334082;C1842036;C2674723;C2750732;CN166718</td>
<td align="left">NA; Epidermal nevus; Congenital giant melanocytic nevus; RAS-associated autoimmune leukoproliferative disorder; Noonan syndrome 6; Rasopathy</td>
</tr>
<tr class="even">
<td align="left">NSD1</td>
<td align="left">64324</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">nuclear receptor binding SET domain protein 1</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">NTHL1</td>
<td align="left">4913</td>
<td align="left">AR</td>
<td align="left">NA</td>
<td align="left">nth like DNA glycosylase 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C4225157</td>
<td align="left">Familial adenomatous polyposis 3</td>
</tr>
<tr class="even">
<td align="left">OGG1</td>
<td align="left">4968</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">8-oxoguanine DNA glycosylase</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">PALB2</td>
<td align="left">79728</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">partner and localizer of BRCA2</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1835817</td>
<td align="left">Fanconi anemia, complementation group N</td>
</tr>
<tr class="even">
<td align="left">PAX5</td>
<td align="left">5079</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">paired box 5</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C3809874</td>
<td align="left">Leukemia, acute lymphoblastic, susceptibility to, 3</td>
</tr>
<tr class="odd">
<td align="left">PDGFRA</td>
<td align="left">5156</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">platelet derived growth factor receptor alpha</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0238198</td>
<td align="left">Gastrointestinal stroma tumor</td>
</tr>
<tr class="even">
<td align="left">PHOX2B</td>
<td align="left">8929</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">paired like homeobox 2b</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C2751682</td>
<td align="left">Neuroblastoma 2</td>
</tr>
<tr class="odd">
<td align="left">PIK3CA</td>
<td align="left">5290</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">PINK1</td>
<td align="left">65018</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">PTEN induced putative kinase 1</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">PMS1</td>
<td align="left">5378</td>
<td align="left">AR</td>
<td align="left">NA</td>
<td align="left">PMS1 homolog 1, mismatch repair system component</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">CN029768</td>
<td align="left">Familial colorectal cancer</td>
</tr>
<tr class="even">
<td align="left">PMS2</td>
<td align="left">5395</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">PMS1 homolog 2, mismatch repair system component</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0009405;C0265325;C1333990</td>
<td align="left">Hereditary nonpolyposis colon cancer; Turcot syndrome; Lynch syndrome</td>
</tr>
<tr class="odd">
<td align="left">POLD1</td>
<td align="left">5424</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">DNA polymerase delta 1, catalytic subunit</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">CN237711</td>
<td align="left">Polymerase proofreading-related adenomatous polyposis</td>
</tr>
<tr class="even">
<td align="left">POLE</td>
<td align="left">5426</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">DNA polymerase epsilon, catalytic subunit</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">CN237711</td>
<td align="left">Polymerase proofreading-related adenomatous polyposis</td>
</tr>
<tr class="odd">
<td align="left">POLH</td>
<td align="left">5429</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">DNA polymerase eta</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C1848410</td>
<td align="left">Xeroderma pigmentosum, variant type</td>
</tr>
<tr class="even">
<td align="left">POLQ</td>
<td align="left">10721</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">DNA polymerase theta</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">POT1</td>
<td align="left">25913</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">protection of telomeres 1</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C4014476</td>
<td align="left">Melanoma, cutaneous malignant, susceptibility to, 10</td>
</tr>
<tr class="even">
<td align="left">PPM1D</td>
<td align="left">8493</td>
<td align="left">Mosaic</td>
<td align="left">LoF</td>
<td align="left">protein phosphatase, Mg2+/Mn2+ dependent 1D</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">PRDM9</td>
<td align="left">56979</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">PR/SET domain 9</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">PRF1</td>
<td align="left">5551</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">perforin 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1863727</td>
<td align="left">Hemophagocytic lymphohistiocytosis, familial, 2</td>
</tr>
<tr class="odd">
<td align="left">PRKAR1A</td>
<td align="left">5573</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">protein kinase cAMP-dependent type I regulatory subunit alpha</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0406810</td>
<td align="left">Carney complex</td>
</tr>
<tr class="even">
<td align="left">PRSS1</td>
<td align="left">5644</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">serine protease 1</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C0238339</td>
<td align="left">Hereditary pancreatitis</td>
</tr>
<tr class="odd">
<td align="left">PTCH1</td>
<td align="left">5727</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">patched 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0004779</td>
<td align="left">Gorlin syndrome</td>
</tr>
<tr class="even">
<td align="left">PTEN</td>
<td align="left">5728</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">phosphatase and tensin homolog</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1959582;CN072330</td>
<td align="left">PTEN hamartoma tumor syndrome; Cowden syndrome 1</td>
</tr>
<tr class="odd">
<td align="left">PTPN11</td>
<td align="left">5781</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">protein tyrosine phosphatase, non-receptor type 11</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C0028326</td>
<td align="left">Noonan syndrome</td>
</tr>
<tr class="even">
<td align="left">PTPN13</td>
<td align="left">5783</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">protein tyrosine phosphatase, non-receptor type 13</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">RAD50</td>
<td align="left">10111</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">RAD50 double strand break repair protein</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">RAD51</td>
<td align="left">5888</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">RAD51 recombinase</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">RAD51B</td>
<td align="left">5890</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">RAD51 paralog B</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">RAD51C</td>
<td align="left">5889</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">RAD51 paralog C</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C3150653</td>
<td align="left">Fanconi anemia, complementation group O</td>
</tr>
<tr class="odd">
<td align="left">RAD51D</td>
<td align="left">5892</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">RAD51 paralog D</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C3280345</td>
<td align="left">Breast-ovarian cancer, familial 4</td>
</tr>
<tr class="even">
<td align="left">RAF1</td>
<td align="left">5894</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">Raf-1 proto-oncogene, serine/threonine kinase</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C1969056;C1969057</td>
<td align="left">LEOPARD syndrome 2; Noonan syndrome 5</td>
</tr>
<tr class="odd">
<td align="left">RB1</td>
<td align="left">5925</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">RB transcriptional corepressor 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0029463;C0035335</td>
<td align="left">Osteosarcoma; Retinoblastoma</td>
</tr>
<tr class="even">
<td align="left">RCC2</td>
<td align="left">55920</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">regulator of chromosome condensation 2</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">RECQL</td>
<td align="left">5965</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">RecQ like helicase</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">RECQL4</td>
<td align="left">9401</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">RecQ like helicase 4</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0032339</td>
<td align="left">Rothmund-Thomson syndrome</td>
</tr>
<tr class="odd">
<td align="left">RET</td>
<td align="left">5979</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">ret proto-oncogene</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0025268;C0025269;C1833921</td>
<td align="left">Multiple endocrine neoplasia, type 2a; Multiple endocrine neoplasia, type 2b; Familial medullary thyroid carcinoma</td>
</tr>
<tr class="even">
<td align="left">RFWD3</td>
<td align="left">55159</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">ring finger and WD repeat domain 3</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">RHBDF2</td>
<td align="left">79651</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">rhomboid 5 homolog 2</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C1835664</td>
<td align="left">Howel-Evans syndrome</td>
</tr>
<tr class="even">
<td align="left">RING1</td>
<td align="left">6015</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">ring finger protein 1</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">RINT1</td>
<td align="left">60561</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">RAD50 interactor 1</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">RMRP</td>
<td align="left">6023</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">RNA component of mitochondrial RNA processing endoribonuclease</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C0220748</td>
<td align="left">Metaphyseal chondrodysplasia, McKusick type</td>
</tr>
<tr class="odd">
<td align="left">RUNX1</td>
<td align="left">861</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">runt related transcription factor 1</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C1832388</td>
<td align="left">Familial platelet disorder with associated myeloid malignancy</td>
</tr>
<tr class="even">
<td align="left">SBDS</td>
<td align="left">51119</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">SBDS, ribosome maturation factor</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0272170</td>
<td align="left">Shwachman syndrome</td>
</tr>
<tr class="odd">
<td align="left">SCG5</td>
<td align="left">6447</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">secretogranin V</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">SDHA</td>
<td align="left">6389</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">succinate dehydrogenase complex flavoprotein subunit A</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C3279992</td>
<td align="left">Paragangliomas 5</td>
</tr>
<tr class="odd">
<td align="left">SDHAF2</td>
<td align="left">54949</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">succinate dehydrogenase complex assembly factor 2</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1708353</td>
<td align="left">Hereditary Paraganglioma-Pheochromocytoma Syndromes</td>
</tr>
<tr class="even">
<td align="left">SDHB</td>
<td align="left">6390</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">succinate dehydrogenase complex iron sulfur subunit B</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1708353;C1847319</td>
<td align="left">Hereditary Paraganglioma-Pheochromocytoma Syndromes; Paraganglioma and gastric stromal sarcoma</td>
</tr>
<tr class="odd">
<td align="left">SDHC</td>
<td align="left">6391</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">succinate dehydrogenase complex subunit C</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1708353;C1847319</td>
<td align="left">Hereditary Paraganglioma-Pheochromocytoma Syndromes; Paraganglioma and gastric stromal sarcoma</td>
</tr>
<tr class="even">
<td align="left">SDHD</td>
<td align="left">6392</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">succinate dehydrogenase complex subunit D</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1708353;C1847319</td>
<td align="left">Hereditary Paraganglioma-Pheochromocytoma Syndromes; Paraganglioma and gastric stromal sarcoma</td>
</tr>
<tr class="odd">
<td align="left">SERPINA1</td>
<td align="left">5265</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">serpin family A member 1</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C0221757</td>
<td align="left">Alpha-1-antitrypsin deficiency</td>
</tr>
<tr class="even">
<td align="left">SETBP1</td>
<td align="left">26040</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">SET binding protein 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86</td>
<td align="left">C0265227</td>
<td align="left">Schinzel-Giedion syndrome</td>
</tr>
<tr class="odd">
<td align="left">SH2B3</td>
<td align="left">10019</td>
<td align="left">AR</td>
<td align="left">NA</td>
<td align="left">SH2B adaptor protein 3</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">SH2D1A</td>
<td align="left">4068</td>
<td align="left">XLR</td>
<td align="left">LoF</td>
<td align="left">SH2 domain containing 1A</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C1868674</td>
<td align="left">Lymphoproliferative syndrome 1, X-linked</td>
</tr>
<tr class="odd">
<td align="left">SHOC2</td>
<td align="left">8036</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">SHOC2, leucine rich repeat scaffold protein</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C1843181</td>
<td align="left">Noonan syndrome-like disorder with loose anagen hair 1</td>
</tr>
<tr class="even">
<td align="left">SLC25A13</td>
<td align="left">10165</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">solute carrier family 25 member 13</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">SLX4</td>
<td align="left">84464</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">SLX4 structure-specific endonuclease subunit</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">SMAD4</td>
<td align="left">4089</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">SMAD family member 4</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0345893</td>
<td align="left">Juvenile polyposis syndrome</td>
</tr>
<tr class="odd">
<td align="left">SMARCA4</td>
<td align="left">6597</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">SWI/SNF related, matrix associated, actin dependent regulator of chromatin, subfamily a, member 4</td>
<td align="left">TCGA_PANCAN_2018,NCGC</td>
<td align="left">C2750074</td>
<td align="left">Rhabdoid tumor predisposition syndrome 2</td>
</tr>
<tr class="even">
<td align="left">SMARCB1</td>
<td align="left">6598</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">SWI/SNF related, matrix associated, actin dependent regulator of chromatin, subfamily b, member 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1836327</td>
<td align="left">Rhabdoid tumor predisposition syndrome 1</td>
</tr>
<tr class="odd">
<td align="left">SMARCE1</td>
<td align="left">6605</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">SWI/SNF related, matrix associated, actin dependent regulator of chromatin, subfamily e, member 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1333989</td>
<td align="left">Meningioma, familial</td>
</tr>
<tr class="even">
<td align="left">SOS1</td>
<td align="left">6654</td>
<td align="left">AD</td>
<td align="left">GoF</td>
<td align="left">SOS Ras/Rac guanine nucleotide exchange factor 1</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C0028326</td>
<td align="left">Noonan syndrome</td>
</tr>
<tr class="odd">
<td align="left">SPINK1</td>
<td align="left">6690</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">serine peptidase inhibitor, Kazal type 1</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">SPOP</td>
<td align="left">8405</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">speckle type BTB/POZ protein</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">SPRED1</td>
<td align="left">161742</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">sprouty related EVH1 domain containing 1</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">SPRTN</td>
<td align="left">83932</td>
<td align="left">AR</td>
<td align="left">NA</td>
<td align="left">SprT-like N-terminal domain</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C4015461</td>
<td align="left">Ruijs-Aalfs syndrome</td>
</tr>
<tr class="odd">
<td align="left">SRY</td>
<td align="left">6736</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">sex determining region Y</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C0206661</td>
<td align="left">Gonadoblastoma</td>
</tr>
<tr class="even">
<td align="left">STAT3</td>
<td align="left">6774</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">signal transducer and activator of transcription 3</td>
<td align="left">TCGA_PANCAN_2018,CGC_86</td>
<td align="left">C2936739</td>
<td align="left">Hyper IgE Syndrome, Autosomal Dominant</td>
</tr>
<tr class="odd">
<td align="left">STK11</td>
<td align="left">6794</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">serine/threonine kinase 11</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0031269</td>
<td align="left">Peutz-Jeghers syndrome</td>
</tr>
<tr class="even">
<td align="left">SUFU</td>
<td align="left">51684</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">SUFU negative regulator of hedgehog signaling</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0025149;C1333989</td>
<td align="left">Medulloblastoma; Meningioma, familial</td>
</tr>
<tr class="odd">
<td align="left">TERF2IP</td>
<td align="left">54386</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">TERF2 interacting protein</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">TERT</td>
<td align="left">7015</td>
<td align="left">AD/AR</td>
<td align="left">LoF</td>
<td align="left">telomerase reverse transcriptase</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1851970</td>
<td align="left">Dyskeratosis congenita autosomal dominant</td>
</tr>
<tr class="odd">
<td align="left">TGFBR1</td>
<td align="left">7046</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">transforming growth factor beta receptor 1</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C0546476</td>
<td align="left">Multiple self healing squamous epithelioma</td>
</tr>
<tr class="even">
<td align="left">TGFBR2</td>
<td align="left">7048</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">transforming growth factor beta receptor 2</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">TMEM127</td>
<td align="left">55654</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">transmembrane protein 127</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0031511</td>
<td align="left">Pheochromocytoma</td>
</tr>
<tr class="even">
<td align="left">TP53</td>
<td align="left">7157</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">tumor protein p53</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1835398</td>
<td align="left">Li-Fraumeni syndrome 1</td>
</tr>
<tr class="odd">
<td align="left">TP63</td>
<td align="left">8626</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">tumor protein p63</td>
<td align="left">CGC_86</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">TRIM37</td>
<td align="left">4591</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">tripartite motif containing 37</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">TSC1</td>
<td align="left">7248</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">TSC complex subunit 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1854465</td>
<td align="left">Tuberous sclerosis 1</td>
</tr>
<tr class="even">
<td align="left">TSC2</td>
<td align="left">7249</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">TSC complex subunit 2</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C1860707</td>
<td align="left">Tuberous sclerosis 2</td>
</tr>
<tr class="odd">
<td align="left">TSHR</td>
<td align="left">7253</td>
<td align="left">AD</td>
<td align="left">NA</td>
<td align="left">thyroid stimulating hormone receptor</td>
<td align="left">TCGA_PANCAN_2018,CGC_86</td>
<td align="left">C1836706</td>
<td align="left">Hyperthyroidism, nonautoimmune</td>
</tr>
<tr class="even">
<td align="left">UROD</td>
<td align="left">7389</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">uroporphyrinogen decarboxylase</td>
<td align="left">TCGA_PANCAN_2018</td>
<td align="left">C0268323</td>
<td align="left">Familial porphyria cutanea tarda</td>
</tr>
<tr class="odd">
<td align="left">VHL</td>
<td align="left">7428</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">von Hippel-Lindau tumor suppressor</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0019562</td>
<td align="left">Von Hippel-Lindau syndrome</td>
</tr>
<tr class="even">
<td align="left">WAS</td>
<td align="left">7454</td>
<td align="left">XLR</td>
<td align="left">LoF</td>
<td align="left">Wiskott-Aldrich syndrome</td>
<td align="left">TCGA_PANCAN_2018,CGC_86</td>
<td align="left">C1839163</td>
<td align="left">Thrombocytopenia, X-linked</td>
</tr>
<tr class="odd">
<td align="left">WRN</td>
<td align="left">7486</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">Werner syndrome RecQ like helicase</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0043119</td>
<td align="left">Werner syndrome</td>
</tr>
<tr class="even">
<td align="left">WT1</td>
<td align="left">7490</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">Wilms tumor 1</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">CN033288</td>
<td align="left">Wilms tumor 1</td>
</tr>
<tr class="odd">
<td align="left">XPA</td>
<td align="left">7507</td>
<td align="left">AR</td>
<td align="left">LoF</td>
<td align="left">XPA, DNA damage recognition and repair factor</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C0268135</td>
<td align="left">Xeroderma pigmentosum, type 1</td>
</tr>
<tr class="even">
<td align="left">XPC</td>
<td align="left">7508</td>
<td align="left">AR</td>
<td align="left">NA</td>
<td align="left">XPC complex subunit, DNA damage recognition and repair factor</td>
<td align="left">TCGA_PANCAN_2018,CGC_86,NCGC</td>
<td align="left">C2752147</td>
<td align="left">Xeroderma pigmentosum, group C</td>
</tr>
<tr class="odd">
<td align="left">XRCC2</td>
<td align="left">7516</td>
<td align="left">AD</td>
<td align="left">LoF</td>
<td align="left">X-ray repair cross complementing 2</td>
<td align="left">NCGC</td>
<td align="left"></td>
<td align="left">NA</td>
</tr>
</tbody>
</table>
