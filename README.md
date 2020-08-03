# thesis-supplementary

*Online supplementary material from the doctoral thesis of Antoine Simon*

The folders `AS` and `SW` contain supplementary material regarding the metatranscriptomic methods and tools described in __Section 2.2__ of the thesis.
The folder `Videos` includes two supplementary videos mentionned in __Section 2.3.3__ of the thesis.


## AS

The dataset consists of twelve libraries treated under two different conditions and using an **poly(A) RNA selection** protocol (i.e., referred to as dataset **AS** in this chapter):

* Cyanomorph (AS1, AS2, AS8, AS9, AS11, AS12)
* Chloromorph (AS3, AS4, AS5, AS6, AS7, AS10)

The R Markdown `AS/scripts.md` file includes 17 scripts used namely for collapsing and filtering [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) transcripts that have been analyzed using [Transdecoder](https://github.com/TransDecoder/TransDecoder) to derive protein (ORF) predictions, analyzed using [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [eXpress](https://pachterlab.github.io/eXpress/) to generate transcriptome-specific abundance estimates. Finally, differential expression between photomorphs was analyzed using [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html). All these steps are succintly described in __Section 2.2__.

Additionally, the `AS` folder contains four supplementary tables:

* `AS/Supplementary_Table_AS1.txt` is a non-redundant matrix with eXpress total read counts and the various fields of the DIAMOND blast tabular file. It corresponds to the output of the [12th script](https://github.com/ant1simon/thesis-supplementary/blob/master/AS/scripts.md#combining-total-read-counts-and-diamond-output).
* `AS/Supplementary_Table_AS2.txt` is the table resulting from the [differential expression analysis](https://github.com/ant1simon/thesis-supplementary/blob/master/AS/scripts.md#performing-differential-expression-analysis). Only **fungal** transcripts with counts per million (CPM) of 20 or greater for at least two samples were included in the analysis. 
* `AS/Supplementary_Table_AS3.txt` contains all statistically-significant **fungal** transcripts with Gene Ontology (GO) terms assigned via the Blast2GO methodology.
* `AS/Supplementary_Table_AS4.txt`includes enriched GO terms (most specific ones only). GO term enrichment was performed on annotations of differentially expressed genes (FDR $<$ 0.05) using Fisher's exact Test.

## SW

The data consists of four libraries treated under two different conditions and using an **rRNA depletion** protocol (i.e., referred to as dataset **SW** in this chapter):

* Cyanomorph (SW1, SW3)
* Chloromorph (SW5, SW9)

The R Markdown `SW/scripts.md` file includes 18 scripts used for collapsing and filtering [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) transcripts that have been analyzed using [Transdecoder](https://github.com/TransDecoder/TransDecoder) to derive protein (ORF) predictions, analyzed using [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [eXpress](https://pachterlab.github.io/eXpress/) to generate transcriptome-specific abundance estimates. Finally, differential expression between photomorphs was analyzed using [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html). All these steps are succintly described in __Section 2.2__.

Additionally, the `SW` folder contains four supplementary tables:

* `SW/Supplementary_Table_SW1.txt` is a non-redundant matrix with eXpress total read counts and the various fields of the DIAMOND blast tabular file. It corresponds to the output of the [12th script](https://github.com/ant1simon/thesis-supplementary/blob/master/SW/scripts.md#combining-total-read-counts-and-diamond-output).
* `SW/Supplementary_Table_SW2.txt` is the table resulting from the [differential expression analysis](https://github.com/ant1simon/thesis-supplementary/blob/master/SW/scripts.md#performing-differential-expression-analysis). Only **bacterial transcripts (except from cyanobacteria)** with counts per million (CPM) of 20 or greater for at least two samples were included in the analysis. 
* `SW/Supplementary_Table_SW3.txt` contains all statistically-significant **bacterial transcripts (except from cyanobacteria)** with Gene Ontology (GO) terms assigned via the Blast2GO methodology.
* `SW/Supplementary_Table_SW4.txt`includes enriched GO terms (most specific ones only). GO term enrichment was performed on annotations of differentially expressed genes (FDR $<$ 0.05) using Fisher's exact Test.

## Videos

* Video S1: Animation of a three-dimensional rendering from confocal z-stacks through a fluorescently hybridized chloromorphic thallus of *Dendriscosticta gelida*, with 6-FAM-labelled Rhizobiales visible in the upper cortex. Green, Rhizobiales; red, algal fluorescence; yellow, *D. gelida*.

* Video S2: Animation of a three-dimensional rendering from confocal z-stacks through a fluorescently hybridized cyanomorphic thallus of *Dendriscosticta gelida*. Blue, DAPI; Green, Rhizobiales; red, algal fluorescence; yellow, *D. gelida*.
