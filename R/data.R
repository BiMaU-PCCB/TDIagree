#' Acute Myeloid Leukaemia agreement data
#'
#' @description
#' Acute myeloid leukaemia (AML) is a type of cancer that starts in the blood-forming cells of the bone marrow.
#' While in adults it is the most common type of leukaemia, it is much rarer in children, accounting for
#' 15-20\% percent of paediatric leukaemia cases, which translates to 8 cases per year for every million children
#' under the age of 15 years.
#'
#' Minimal residual disease (MRD) is the percentage of cancer cells that remain in a person
#' either during or after treatment when the patient is in remission (no symptoms or signs of disease). MRD aids
#' in identifying high-risk patients so therapy can be intensified in them while deintensification of
#' therapy can prevent long-term sequelae of chemotherapy in low-risk category patients.
#'
#' MRD describes disease that can be detected using techniques other than traditional morphology, including
#' molecular methods such as polymerase chain reaction (PCR) and immunological methods such as flow cytometry
#' (FCM) (Chattaerjee \emph{et al.}, 2016).
#'
#' This dataset is adapted from the \emph{Childhood Leukemia: Overcoming distance between South America
#' and Europe Regions} (CLOSER) project, whose goal was to decrease the gap between Europe and Latin America in terms
#' of the diagnosis, monitoring, survival, and quality of life of patients with childhood leukaemia and their caregivers.
#' See \strong{Source} for further information on the project. The dataset contains data from 116 paediatric patients
#' diagnosed with AML, in which the MRD was measured twice after treatment initiation by the methods PCR and FCM.
#'
#' @format
#'
#' A data frame in long format with the following columns:
#'
#' \tabular{ll}{
#' id: \tab Patient identifier \cr
#' met: \tab Method to quantify MRD (PCR or FCM) \cr
#' rep: \tab Replicate (1 = first, 2 = second) \cr
#' mrd: \tab MRD (\%) \cr
#' }
#'
#'
#' @references Chatterjee, T., Mallhi, R. S., & Venkatesan, S. (2016). Minimal residual disease detection using flow cytometry: Applications in acute leukemia. Medical Journal Armed Forces India, 72(2), 152-156.
#'
#' @source \url{https://closerleukemia.eu/}
"AMLad"
