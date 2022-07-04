# Read me file

Manuscript: “20-Year Benefit from Adjuvant Goserelin and Tamoxifen in Premenopausal Breast Cancer Patients in a Controlled Randomized Clinical Trial”</br>

<b>Kaplan_Meier_Figure2_Figure3.R</b></br>
Description: Performs Kaplan-Meier analysis on:</br>
•	Goserelin vs Control</br>
•	Tamoxifen vs Control</br>
•	Goserelin+Tamoxifen vs Control</br>
Unstratified (all patients) and stratified analyses based on patients’ genomic risk (70-gene signature low or high risk)</br>
Input: STO5_trial.RData</br>
Output: Figure 2, Figure 3</br>

<b>Cox_analyses_Figure2_Figure3.R</b></br>
Description: Performs multivariable Cox proportional hazard regression on:</br>
•	Goserelin vs Control</br>
•	Tamoxifen vs Control</br>
•	Goserelin+Tamoxifen vs Control</br>
Unstratified (all patients) and stratified analyses based on patients’ genomic risk (70-gene signature low or high risk)</br>
Input: STO5_trial.RData</br>
Output: Hazard Ratios and 95% CIs for Figure 2 and Figure 3</br>

<b>Cox_analyses_interaction_Table2.R</b></br>
Description: Performs multivariable Cox proportional hazard regression to analyze:</br>
1.	the effect of Goserelin in patients treated with and without Tamoxifen:</br>
a.	Goserelin+Tamoxifen vs Tamoxifen</br>
b.	Goserelin vs Control</br>
2.	the effect of Tamoxifen in patients treated with and without Goserelin:</br>
a.	Goserelin+Tamoxifen vs Goserelin</br>
b.	Tamoxifen vs Control</br>
Also performs interaction test between Goserelin and Tamoxifen by including a product term in the Cox model.</br>
Input: STO5_trial.RData</br>
Output: Table 2</br>

<b>Time_varying_analysis_SFigure1_Table3.R</b></br>
Description: Performs multivariable parametric flexible modelling on:</br>
•	Tamoxifen vs Control in genomic low-risk patients</br>
•	Goserelin vs Control in genomic high-risk patients.</br>
Input: STO5_trial.RData</br>
Output: Supplementary Figure S1 and Table 3 with estimated hazard ratios at years 5, 10, 15 and 20.</br>
