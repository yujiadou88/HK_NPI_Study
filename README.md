# Hong Kong NPI study

## 1. Study overview {#study-overview}

There is a serious deficit in the evidence base on whether non-pharmaceutical interventions can reduce the spread of influenza. We implemented a study of whether face masks and hand hygiene can reduce influenza transmission among Hong Kong household members.

In 2007 (pilot) and 2008 (main study) we recruited subjects presenting to outpatient clinics (in both the private and public sectors across Hong Kong) with influenza-like-illness of \<48 hours duration. After influenza was confirmed in an index case by the QuickVue Influenza A+B rapid test, or in some clinics based on a clinical definition alone, the household of the index subject was randomized to 1) control or 2) hand hygiene or 3) hand hygiene plus surgical face masks (the 3rd arm was surgical face masks in the pilot study). We aimed to implement the interventions at an initial home visit within 36 hours of recruitment, and evaluate subsequent infections by self-reported daily symptom diaries and home visits after 3 and 6 days (and additionally after 9 days in the pilot study). Nose and throat swabs were collected from index subjects and all household contacts at each home visit and tested by viral culture and or RT-PCR. The primary outcome measure was laboratory-confirmed influenza infection in a household contact by viral culture (pilot) or RT-PCR (main); the secondary outcome was clinically diagnosed influenza by self-reported symptoms. We evaluated adherence by self-report and by measuring the number of surgical masks used and weighing the amount of soap/hand rub used.

Full details of our study design are available in our [study protocol](https://doi.org/10.1371/journal.pone.0002101) (the link is to version GML001.5, 10-Dec-2007).

We are often asked why we chose to recruit index cases presenting with ILI and observe their households for a week (called a 'case-ascertained' design), rather than starting with a cohort of uninfected households and following them over an influenza season. Much greater resources might be required for the latter study, given the low annual attack rate of influenza.

A weakness of our chosen study design is the delay between index case symptom onset (or more precisely the onset of index case infectiousness) and application of the intervention, and any intervention effects we observe would likely be attenuated due to this delay. We are actively working on models which can estimate intervention effects allowing for delays. Key changes made:

## 2. Raw data

The latest version of the NPI study year 1 (pilot) data (02-2007 to 09-2007) are available to download as a zip file here:

-   [HongKongNPIpilot.zip](data/HongKongNPIpilot.zip).

This version of the dataset covers the transmission of influenza in households and adherence to interventions, antiviral treatment, quantitative viral loads, and data from recruiting clinics including presenting symptoms and rapid test performance in 2007.

The latest version of the NPI study year 2 data (01-2008 to 09-2008) are available to download as a zip file here:

-   [HongKongNPIstudy.zip](data/HongKongNPIstudy.zip).

This version of the dataset covers the transmission of influenza in households and adherence to interventions, antiviral treatment, quantitative viral loads, and data from recruiting clinics including presenting symptoms and rapid test performance in 2008.

The latest version of the NPI study year 3 data (01-2009 to 06-2009) are available to download as a zip file here:

-   [HongKongNPIstudy2009V1.zip](data/HongKongNPIstudy2009V1.zip).

This version of the dataset covers the transmission of influenza in households and subject demographics, antiviral treatment, quantitative viral loads, and data from recruiting clinics including presenting symptoms and rapid test performance in 2009.

We provide our data under the [Open Data Commons Public Domain Dedication and License](http://www.opendatacommons.org/odc-public-domain-dedication-and-licence/), which is a version of open access for data. Under this licence we reserve no rights: there are no restrictions on use of our data, and no requirement to cite our work or this website. However we would anticipate that for academic purposes the standard practice of referencing sources would apply. We would like to hear from researchers who are using our data and we would be keen to work together on analyses.

## 3. Pilot findings

The primary findings from the 2007 pilot study were published by [Cowling et al. (2008, PLoS ONE)](http://dx.doi.org/10.1371/journal.pone.0002101 "full text of Cowling et al., 2008, PLoS ONE"). In brief, the pilot study achieved its aim of confirming the feasibility of the main study design, and provided some information about the secondary attack ratios. We found no evidence of intervention effects although with only 122 households our pilot study was underpowered to detect even fairly large differences between arms, and adherence to the hand hygiene and face mask interventions was variable.

All results can be reproduced using these R scripts:

-   [Table 1](NPIpilot_scripts/Table_1.r).

-   [Table 2](NPIpilot_scripts/Table_2.r).

-   [Table 3](NPIpilot_scripts/Table_3.r) *(requires R package 'gee')*.

-   [Figure 1](NPIpilot_scripts/Figure_1.r).

-   [Figure 2](NPIpilot_scripts/Figure_2.r).

-   [Table S1](NPIpilot_scripts/Table_S1.r) *(requires R package 'ROCR')*.

-   [Figure S1](NPIpilot_scripts/Figure_S1.r).

-   [Text results](NPIpilot_scripts/Text_results.r).

## 4. Serial interval

We investigated the clinical onset serial intervals in our pilot study, adjusting for our 'case-ascertained' study design, and our results are published in [Cowling et al. (2009, Epidemiology)](http://dx.doi.org/10.1097/EDE.0b013e31819d1092). We used symptom onset times from 14 pairs of infector/infectee in the NPI pilot study to estimate the clinical serial interval of human influenza in households to have mean 3.6 days (95% confidence interval: 2.9, 4.3 days), with standard deviation 1.6 days. This is slightly longer than some previous estimates.

Our finding of a mean serial interval of 3.6 days suggests that the average time from symptom onset in the index case to secondary infection in the household setting may be around 2 days, assuming that the average time from secondary infection to secondary onset (the incubation period) is around 1.5 days ([Moser et al., AJE, 1979](http://www.ncbi.nlm.nih.gov/pubmed/463858)). In the context of our NPI study where index cases were recruited within 48 hours of symptom onset, these findings suggest that our NPI study design may lead to an underestimate of the true effectiveness of the interventions because some infections may have occurred prior to recruitment and intervention. Nevertheless, provided that interventions can be applied soon after symptom onset, it is likely that in our NPI study we would be able to observe attenuated efficacies, and this should be taken into consideration when interpreting the results.

Results described in Cowling et al. (2009, Epidemiology) are reproduced in the following scripts which can be run in [R](http://www.r-project.org):

-   [functions](Serial_scripts/functions.r) to calculate the non-parametric and parametric estimates of the serial interval (used in some of the other scripts here).

-   [Figure 1](Serial_scripts/Figure_1.r).

-   [Figure 2](Serial_scripts/Figure_2.r).

-   [Appendix Table 1](Serial_scripts/Appendix_Table_1.r) which produces the dataset [serial_appendix_table_1.csv](Serial_scripts/serial_appendix_table_1.csv) used in some of the other scripts here.

-   [Appendix Table 2](Serial_scripts/Appendix_Table_2.r).

-   [Appendix Figure 1](Serial_scripts/Appendix_Figure_1.r).

-   Appendix Table 3 and Appendix Figure 2 *(downloadable zip file with WinBUGS scripts - these require the R package R2WinBUGS and a bit of effort setting up WinBUGS to run in R - instructions available [here](http://www.stat.columbia.edu/~gelman/bugsR/))*.

## 5. Rapid test sensitivity

We examined the performance of the QuickVue Influenza A+B rapid test, which we used in our study to screen subjects. The results are published in [Cheng et al. (2009, Diagnostic Microbiology and Infectious Disease)](http://dx.doi.org/10.1016/j.diagmicrobio.2009.05.003). We studied the sensitivity and specificity of the rapid test versus viral culture as the gold standard. We found moderate overall sensitivity of around 68% and high specificity of around 96%.

We also tested a subset of specimens by quantitative RT-PCR and found that test sensitivity increases with viral load. We also found that viral loads were higher in children and sooner after illness onset.

Results described in Cheng et al. (2009) are reproduced in the following scripts which can be run in [R](http://www.r-project.org):

-   [functions](QuickVue_scripts/qv_functions.r) to calculate sensitivity, specificity,and 95% confidence intervel (used in some of the other scripts here).

-   [scripts](QuickVue_scripts/add_groups.r) to reformat the raw data (used in some of the other scripts here).

-   [Table 1](QuickVue_scripts/Table_1.r).

-   [Table 2](QuickVue_scripts/Table_2.r).

-   [Figure 2](QuickVue_scripts/Figure_2.r).

-   [Appendix Table A](QuickVue_scripts/Appendix_Table_A.r).

-   [Appendix Figure 1](QuickVue_scripts/Appendix_Figure_1.r).

-   [Other results described in the manuscript text](QuickVue_scripts/Text_results.r).

## 6. Main results

The primary findings from our main study in 2008 have been published in [Cowling et al. (2009, Annals of Internal Medicine)](http://www.annals.org/cgi/content/full/151/7/437). Our results suggested that hand hygiene and facemasks could prevent a substantial proportion of household transmission of influenza virus when implemented within 36 hours of illness onset in the index case.

Results described in Cowling et al. (2009) are reproduced in the following scripts which can be run in [R](http://www.r-project.org):

-   [scripts](NPImain_scripts/Analyzed_hh.r) to mark the analyzed households (used in some of the other scripts here).

-   [Table 2](NPImain_scripts/Table_2.r).

-   [Table 3](NPImain_scripts/Table_3.r).

-   [Table 4](NPImain_scripts/Table_4.r) *(requires the R package 'gee')*.

-   [Table 5](NPImain_scripts/Table_5.r) *(requires the R package 'gee')*.

-   [Table 6](NPImain_scripts/Table_6.r).

-   [Figure 1](NPImain_scripts/Figure_1.r).

-   [Appendix Table 1](NPImain_scripts/Appendix_table_1.r).

-   [Appendix Table 2](NPImain_scripts/Appendix_table_2.r).

-   [Appendix Table 3](NPImain_scripts/Appendix_table_3.r).

-   [Appendix Table 4](NPImain_scripts/Appendix_table_4.r).

-   [Appendix Table 5](NPImain_scripts/Appendix_table_5.r).

-   [Appendix Table 6](NPImain_scripts/Appendix_table_6.r) *(requires the R package 'gee')*.

-   [Appendix Table 7](NPImain_scripts/Appendix_table_7.r) *(requires the R package 'gee')*.

-   [Appendix Table 8](NPImain_scripts/Appendix_table_8.r) *(requires the R packages 'Hmisc' and 'gee')*.

-   [Appendix Table 9](NPImain_scripts/Appendix_table_9.r) *(requires the R packages 'Hmisc' and 'gee')*.

-   [Appendix Table 10](NPImain_scripts/Appendix_table_10.r).

-   [Appendix Figure 1](NPImain_scripts/Appendix_figure_1.r).

-   [Appendix Figure 2](NPImain_scripts/Appendix_figure_2.r).

in [Cowling et al. (2009, Annals of Internal Medicine)](http://www.annals.org/cgi/content/full/151/7/437).

## 7. Direct and indirect benefits of oseltamivir treatment

Almost 25% of index cases in our pilot and main studies were prescribed oseltamivir treatment by their treating physician. We studied the effectiveness of oseltamivir treatment in reducing duration of illness and viral shedding, as well as the indirect benefits in reducing transmission to household members. The findings from our study are published in [Ng et al. (2010, CID)](https://doi.org/10.1086/650458 "full text of Ng et al., 2010, CID"). We found that oseltamivir treatment given within 24 hours reduced the duration of illness by around half (statistically significant), and the risk of transmission to household contacts by around half (not statistically significant).

Results described in Ng et al. (2010) are reproduced in the following scripts which can be run in [R](http://www.r-project.org "R statistical software homepage"):

-   [scripts](Oseltamivir_scripts/dataframe.r) to reformat the raw data (used in some of the other scripts here).

-   [Table 1](Oseltamivir_scripts/Table_1.r).

-   [Table 2](Oseltamivir_scripts/Table_2.r) *(requires the R packages 'chron' and 'survival')*.

-   [Table 3](Oseltamivir_scripts/Table_3.r) *(requires the R packages 'chron' and 'survival')*.

-   [Table 4](Oseltamivir_scripts/Table_4.r) *(requires the R packages 'chron' and 'gee')*.

-   [Figure 1](Oseltamivir_scripts/Figure_1.r) *(requires the R package 'survival')*.

-   [Figure 2](Oseltamivir_scripts/Figure_2.r) *(requires the R package 'chron')*.

-   [Appendix Table A1](Oseltamivir_scripts/Appendix_Table_A1.r).

-   [Appendix Table A2](Oseltamivir_scripts/Appendix_Table_A2.r) *(requires the R packages 'chron' and 'gee')*.

## 8. Viral shedding and clinical illness

In acute influenza virus infections most viral shedding occurs within a few days of illness onset, and the degree of viral shedding correlates with symptoms and tympanic temperature.

Results described in [Lau et al. (2010, JID)](http://dx.doi.org/10.1086/652241 "full text of Lau et al., 2010, JID") are reproduced in the following scripts which can be run in [R](http://www.r-project.org "R statistical software homepage"):

-   [scripts](vshed_scripts/JID_dataframe.r) to reformat the raw data (used in some of the other scripts here).

-   [Table 1](vshed_scripts/Table_1.r).

-   [Figure 1](vshed_scripts/Figure_1.r).

-   [Figure 2](vshed_scripts/Figure_2.r).

-   [Figure 3](vshed_scripts/WinBUGS.zip) *(downloadable zip file with WinBUGS scripts - these require the R package R2WinBUGS and a bit of effort setting up WinBUGS to run in R - instructions available [here](http://www.stat.columbia.edu/~gelman/bugsR/))*.

-   [Appendix Table 1](vshed_scripts/Appendix_Table_1.r).

-   [Appendix Figure 1](vshed_scripts/Appendix_Figure_1.r).

-   [Appendix Figure 2](vshed_scripts/Appendix_Figure_2.r).

## 9. Heterogeneity in viral shedding among seasonal influenza A virus infections

Between January 2008 and August 2009, 637 patients had influenza A virus detected by a rapid test. [Lau et al. (2013, JID)](http://jid.oxfordjournals.org/content/early/2013/02/20/infdis.jit034.long) studies these subjects and found that viral shedding in more heterogeneous in children than adults. The top 20% most infectious subjects were estimated to be responsible for 80%-90% of the total infectiousness.

Results described in the manuscript are reproduced in the following scripts which can be run in [R](http://www.r-project.org "R statistical software homepage"):

## 10. Modes of transmission of influenza A virus in households

Influenza A viruses are believed to spread between humans through contact, large respiratory droplets and small particle droplet nuclei (aerosols), but the relative importance of each of these modes of transmission is unclear. Cowling and colleagues modelled data of 822 subjects recruited from NPI study 2008-2009 to estimate that aerosol transmission accounts for approximately half of all influenza transmission events.

Results described in the manuscript are reproduced in the following scripts which can be run in [R](http://www.r-project.org "R statistical software homepage"):

-   [scripts](NPI_NatureComm_mot_scripts/dataframe.r) to reformat the raw data (used in some of the other scripts here).

-   [Table 1](NPI_NatureComm_mot_scripts/Table_1.r).

-   [Figure 1](NPI_NatureComm_mot_scripts/Figure_1.r).

-   [Figure 2](NPI_NatureComm_mot_scripts/Figure_2.r).

\*Scripts for results presented in Supplementary file will be uploaded soon.

## 10B. Modes of transmission of influenza B virus in household

Continued with the above study on influenza A, we further explore the modes of transmission of influenza B virus based on the same household settings. We found aerosol transmission may be an important mode of spread of influenza B virus. The point estimates of aerosol transmission were slightly lower for influenza B virus compared to previously published estimates for influenza A virus in both Hong Kong and Bangkok.

Results described in the manuscript are reproduced in the following scripts which can be run in [R](http://www.r-project.org "R statistical software homepage"):

-   [scripts](NPI_PLoSOne_motB_scripts/dataframe.r) to reformat the raw data (used in some of the other scripts here).

-   [Table 2](NPI_PLoSOne_motB_scripts/Table_2.r).

-   [Table 4](NPI_PLoSOne_motB_scripts/Table_4.r).

-   [Figure 1](NPI_PLoSOne_motB_scripts/Figure_1.r).

-   [Figure 2](NPI_PLoSOne_motB_scripts/Figure_2.r).

## 11. Further work

We have a series of further analyses underway using the data from our study. More details will follow later.

## Authors and investigators

The principal investigators of this study are Gabriel Leung and [Ben Cowling](https://sph.hku.hk/en/Biography/Cowling-Benjamin-John). The full list of investigators is given in our study protocol (linked in [section 1](#study-overview)). The data were uploaded by [Ben Cowling](https://sph.hku.hk/en/Biography/Cowling-Benjamin-John), and the scripts were written by Vicky Fang and [Ben Cowling](https://sph.hku.hk/en/Biography/Cowling-Benjamin-John).

## A comment on reproducible research

We fully support the [increasing calls](http://dx.doi.org/10.1097/EDE.0b013e318196784a) from the academic community for [epidemiologic analyses to be reproducible](http://dx.doi.org/10.1093/aje/kwj093 "Peng et al., 2006, AJE"), and [raw data from randomized controlled trials to be published](http://dx.doi.org/10.1186/1745-6215-10-17 "Hrynaszkiewicz and Altman, 2009, Trials"), as a part of the wider scientific endeavour to replicate results. Another example of this recommendation is in the [Good Practice Guide for Quantitative Veterinary Epidemiology](http://www.qve-goodpracticeguide.org.uk/guide#TOC-2.4-Inputs). Here we have published the raw *anonymised* data from our HK NPI pilot study, and will later release the data from our main study. We have also published scripts which allow the analyses in our published papers to be reproduced.

Thousands of local people have given their time, and their families', as part of their participation in our studies, all in the expectation that our research studies will add to medical and scientific knowledge. Participants should also expect that we will make the best possible use of the information that we have collected about them. It would be difficult to argue that facilitating best use of the data by the research community need not involve releasing raw data.

Publication of anonymised raw data has been approved by our local IRB and funding sources, and participants were advised that anonymised data would be published during the informed consent process. We anticipate that release of the raw data will:

-   Promote reproducibility of our results.

-   Allow other investigators to conduct their own analyses on our data.

-   Allow other investigators to compare our data with theirs, for example to explore similarities and differences between research findings.

-   Allow other investigators to prepare and plan for their own studies.

## Publications

1.  Cowling BJ, Fung RO, Cheng CK, Fang VJ, Chan KH, Seto WH, Yung R, Chiu B, Lee P, Uyeki TM, Houck PM, Peiris JS, Leung GM. Preliminary findings of a randomized trial of non-pharmaceutical interventions to prevent influenza transmission in households. *PLoS ONE*, 2008; **3**(5):e2101. [[link]](http://dx.doi.org/10.1371/journal.pone.0002101) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/18461182).

2.  Cowling BJ, Fang VJ, Riley S, Peiris JS, Leung GM. Estimation of the serial interval of influenza. *Epidemiology*, 2009; **20**(3):344-7. [[link]](http://dx.doi.org/10.1097/EDE.0b013e31819d1092) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/19279492).

3.  Cheng KYC, Cowling BJ, Chan KH, Fang VJ, Seto WH, Yung R, Uyeki TM, Houck PM, Peiris JSM, Leung GM. Factors affecting QuickVue Influenza A+B rapid test performance in the community setting. *Diagnostic Microbiology and Infectious Disease*, 2009; **65**(1):35-41. [[link]](http://dx.doi.org/10.1016/j.diagmicrobio.2009.05.003) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/19679233).

4.  Cowling BJ, Chan KH, Fang VJ, Cheng CKY, Fung ROP, Wai W, Sin J, Seto WH, Yung R, Chu DWS, Chiu BCF, Lee PWY, Chiu MC, Lee HC, Uyeki TM, Houck PM, Peiris JSM, Leung GM. Facemasks and hand hygiene to prevent influenza transmission in households: a randomized trial. *Annals of Internal Medicine*, 2009; **151**(7):437-46. [[link]](http://www.annals.org/cgi/content/full/151/7/437) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/19652172).

5.  Ng S, Cowling BJ, Fang VJ, Chan KH, Ip DKM, Cheng CKY, Uyeki TM, Houck PM, Peiris JSM, Leung GM. Effects of oseltamivir treatment on duration of clinical illness and viral shedding, and household transmission of influenza virus. *Clinical Infectious Diseases*, 2010; **50**(5):707-14. [[link]](https://doi.org/10.1086/650458) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/20121573).

6.  Lau LLH, Cowling BJ, Fang VJ, Chan KH, Lau EHY, Lipsitch M, Cheng CKY, Mouck PM, Uyeki TM, Peiris JSM, Leung GM. Viral shedding and clinical illness in naturally acquired influenza virus infections. *Journal of Infectious Diseases*, 2010; 201:1509-16. [[link]](http://dx.doi.org/10.1086/652241) [[PubMed]](http://www.ncbi.nlm.nih.gov/pubmed/20377412).

7.  Lau LLH, Ip DKM, Nishiura H, Fang VJ, Chan KH, Peiris JSM, Leung GM, Cowling BJ. Heterogeneity in viral shedding among individuals with medically attended influenza A virus infection. *Journal of Infectious Diseases*, 2013; **207**(8):1281–1285. [[link]](https://academic.oup.com/jid/article/207/8/1281/891314) [[Pubmed]](https://pubmed.ncbi.nlm.nih.gov/23382573/).

## Acknowledgements

The NPI study was funded by the US Centers for Disease Control and Prevention (grant no. 1 U01 CI000439-01). The work presented here has also received funding from the Research Fund for the Control of Infectious Disease, Food and Health Bureau, Government of the Hong Kong SAR (grant nos. 08070632, 08070562 and HKU-AA-23), the Area of Excellence Scheme of the Hong Kong University Grants Committee (grant no. AoE/M-12/06), the US National Institutes of Health cooperative agreement (grant no. 5 U01 GM076497), the Harvard Center for Communicable Disease Dynamcs from the US National Institutes of Health Models of Infectious Disease Agent Study program (grant no. 1 U54 GM088558), and the HKU University Research Council Strategic Research Theme of Public Health.

------------------------------------------------------------------------

[![Creative Commons License](https://i.creativecommons.org/l/by/3.0/80x15.png)](http://creativecommons.org/licenses/by/3.0/) This work is licensed under a [Creative Commons Attribution 3.0 Unported License](http://creativecommons.org/licenses/by/3.0/). It was written by [Ben Cowling](http://www.hku.hk/cmd/staff/bio/cowling.htm) (email).
