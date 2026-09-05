# MS Student in Artificial Intelligence
#### Technical Skills: Python, R, C, SQL


## Education
- MS in Artificial Intelligence, AI and Biomedical concentration | Columbia University   
  *Sep 2026 - Dec 2027*
- BS in Data Science | Ewha Womans University  
  *Mar 2023 - Jun 2026*
- Exchange Student | Western Washington University  
  *Jan 2025 - Jun 2025*
- Summer Session | University of California, Los Angeles  
  *Jun 2023 - Aug 2023*


## Publications, Presentations, Posters
- Du, H., **Cheong, J.**, Zhu, G., Huang, X., Zheng, J., Zhong, T., Hou, C., & Li, S. (Accepted).  
  *What Draws the Upvotes? Examining Relationships Among Performance, Interpersonal Skills, and Learning Behaviors in Online Social Annotation Environments.*  
  Poster session presented at the **American Educational Research Association (AERA)** 2026 Annual Meeting, Los Angeles, CA.  
  📄 [Click to View PDF](assets/pdf/AERA2026.pdf)
  
- Du, H., **Cheong, J.**, Zhu, G., Huang, X., Li, S., Zheng, J., Tan, J., Zhong, T., Hou, C., & Liu, T. (Accepted).  
  *Identifying Student Clusters Based on Interpersonal Skills and Examining the Correlations.*  
  Presentation as a poster at the **19th International Conference on Computer-Supported Collaborative Learning (CSCL 2026)**, International Society of the Learning Sciences.  
  📄 [Click to View PDF](assets/pdf/ISLS2026.pdf)      
 


## Research Experiences
**Undergraduate Research Intern with Professor Hanxiang Du**  
*Apr 2025 - Sep 2025*
- Submitted poster to *AERA 2026 Annual Meeting* 
- Submitted short paper to *CSCL 2026 International Conference*
- Conducted educational data mining, hierarchical clustering, cosine similarity analysis, CRQA
- Studied about NLP and used NLTK for text preprocessing

**Undergraduate Research Intern with Professor Jenny Jyoung Lee**  
*Sep 2024 - Dec 2024*
- Conducted a literature review on air pollution disparities in relation to mobility
- Calculated segregation measures incorporating spatial patterns in census tract using U.S. Census data

## Work Experiences
**Tutor for Introduction to Linear Algebra**  
*Ewha Womans University, Seoul, South Korea (Sep 2025 - Dec 2023)*
- Provided weekly tutoring in English to support students’ understanding of course material  

**Teaching Assistant (Part-Time)**  
*Star Math Academy, Gyunggi, South Korea (Sep 2023 - Dec 2023)*
- Assisted students with understanding and solving math questions and graded daily quiz and homework

**Teaching Assistant (Part-Time)**  
*English Vine Academy, Gyunggi, South Korea (Dec 2021 - Feb 2022)*
- Supported students by delivering further explanation and graded vocabulary quiz and homework

## Projects
### Seat Occupancy Prediction–Based Tiered Pricing Policy for Performance Venues
*Data Science Capstone Design 2, Department of Data Science, Ewha Womans University (Spring 2026) — Advisor: Prof. Sangyeop Kim*  
South Korea's performing arts market has grown steadily, yet individual productions increasingly struggle with profitability, partly because venues still rely on a rigid, 3–5 tier fixed pricing structure. Since tickets are time-limited, capacity-constrained goods, we treated seat pricing as a revenue-management problem and asked whether seat-level demand could be predicted and used to design a more responsive pricing policy. 
Using 25 performances (May 5–24) of the musical *Lempicka* as a case study, we built a seat × performance level dataset by crawling seat-sale snapshots from Interpark NOL, audience reviews from Seeya!, and cast Instagram follower counts. Review text was scored along four dimensions (view, comfort, lighting, sound) and combined with seat location, price tier, cast, and schedule features. Seven classifiers were compared, and CatBoost was selected as the final model (Test F1 = 0.84, ROC-AUC = 0.87); stage distance, performance date, row number, and distance from center emerged as the strongest predictors of whether a seat would sell. 
Using the model's predicted sale probabilities, we simulated an expected-revenue-maximizing pricing policy that produced 7–8 finer price tiers per performance instead of the venue's fixed grid. Across all 25 performances, the proposed policy raised average seat occupancy from 58.8% to 70.6% (+11.75pp) and total revenue by 4.2% (about ₩79 million) relative to actual pricing — suggesting that seat-level, data-driven pricing can improve both occupancy and revenue simultaneously. 
🧑‍🏫 [Click to View Slides (Korean Original)](assets/projects/Ewha2026Spring/Capstone/FinalSlides)  
📄 [Click to View PDF (Korean Original)](assets/projects/Ewha2026Spring/Capstone/FinalReport)
🖼️ [Click to View Poster (Korean Original)](assets/projects/Ewha2026Spring/Capstone/FinalPoster)

### A SigLIP2-Based Adapter Embedding Channel for Long-Tail Personal Gallery Search
*Term project for Deep Learning, Ewha Womans University (Spring 2026)* 
CLIP-style vision-language models retrieve well for object-centric queries, but real users searching their own photo galleries often type long-tail, descriptive queries — "an exciting party photo," "a sentimental night walk" — that general-purpose VLMs like SigLIP2 struggle to align with the right image. Since no large-scale query-answer dataset exists for these personal, mood- or scene-driven queries, we used public image-caption datasets (Flickr30k, COCO) as a proxy to first improve object- and scene-level image-text alignment, freezing the pretrained SigLIP2 encoder and training a lightweight adapter on top of it. 
We explored two directions. A reranker that reorders SigLIP2's top retrieved candidates only helped in narrow settings (Top-5/Top-10, when the model was already uncertain) and was not adopted. A residual MLP adapter (2.66M params), trained with an InfoNCE loss and hard-negative mining (remining the top-K hardest same-batch negatives every 5 epochs) on 620K image-caption pairs from Flickr30k and COCO, proved more effective — it outperformed three alternative adapter architectures (symmetric, bottleneck, attention) in stability and accuracy. 
With tuned batch size/temperature and full 5-caption COCO training, the final model raised Recall@1 from 87.80% to 91.00% on Flickr30k (+3.20pp) and from 49.48% to 55.26% on COCO 5K (+5.78pp) over the frozen SigLIP2 baseline. On a real personal gallery of 964 photos with 7 descriptive queries, it reached 57.1% R@1 and 71.4% R@5 — succeeding on visually distinctive scenes (a wedding aisle, bowls of ramen) but struggling when the gallery held many visually similar photos (several night cityscapes), pointing toward future work that pairs this adapter with metadata- and caption-based retrieval channels for full long-tail search. 
🧑‍🏫 [Click to View Slides (Korean Original)](assets/projects/Ewha2026Spring/DeepLearning/FinalSlides.pdf)

### A VLM-Based Natural Language Long-Tail Fashion Attribute Retrieval System
*Course project for Data Engineering, Ewha Womans University (Spring 2026)* 
Fashion e-commerce search typically struggles with compound natural-language queries combining multiple attributes at once (e.g., "a loose lavender knit long-sleeve") or filtering bottoms by waist rise and fit. We built a long-tail fashion attribute search system to close this gap, using 100,000 polygon-masked images from the K-Fashion dataset (AI Hub).  
Each image was processed through three channels: a frozen marqo-fashionSigLIP visual embedding indexed in Qdrant; attribute classifier heads (pattern, trim, length, waist rise) trained on a shared backbone via iterative active learning across five annotation rounds, reaching 0.73–0.82 accuracy/F1; and VLM-generated captions (Gemini 2.5 Flash Lite) capturing fine-grained details and mood/TPO. A Korean query is translated and corrected against a fashion-term dictionary, then run through visual search with Qdrant filtering, boosted by caption keyword matches, and reranked via a Korean sentence-embedding model (ko-sroberta) blending visual and semantic similarity. 
The final pipeline reached Precision@10 of 0.847 across easy/medium/hard query sets — notably 0.78 vs. 0.46 for a visual-only baseline on medium-difficulty compound queries — versus 0.693 for a text-filtering-only variant. We also built an end-to-end Prefect pipeline routing new images through captioning, classification, and embedding, with a human-in-the-loop Gradio annotation step. 
🧑‍🏫 [Click to View Slides](assets/projects/Ewha2026Spring/DataEngineering/Prj2FinalSlides.PDF)

### Reactions to AI: Korea vs. Western Media and Communities
*Course project for Data Engineering, Ewha Womans University (Spring 2026)* 
Motivated by AI's rapid growth and cultural/linguistic differences in how it's discussed, this project built an end-to-end pipeline comparing Korean and English-speaking reactions to AI across news media, online communities, and professional fields — covering crawling, preprocessing, database storage, text/sentiment analysis, and visualization. 
I was responsible for the crawling stage with a teammate, building source-specific crawlers: requests + BeautifulSoup for static HTML (Nate Pann), Playwright for JS-heavy/login-gated sites (Everytime, Blind, DC Inside, Naver News, Reddit), and official APIs where available (YouTube Data API v3, Guardian API, Arctic Shift for large-scale occupation-based Reddit sampling across 9 professional subreddits). We solved practical scraping challenges including multi-tier selector fallback for inconsistent news HTML, automated login/infinite-scroll handling, exponential backoff for Reddit rate limits, checkpoint/resume for daily API quotas, and batched incremental saves. We collected 19,506 posts across eight core Korean/English sources plus 95,382 occupation-based Reddit posts, stored in MongoDB Atlas after deduplication and tagging. Analysis used BERTopic for topic extraction and multilingual BERT for sentiment classification, visualized with Plotly. 
Key findings: negative sentiment toward AI rose over time in both regions' news (Korea: 44.7%→66.7%; Western: 54.7%→70.9%), with Korean discourse shifting from growth-engine framing toward geopolitical risk (DeepSeek, China, U.S. policy), and Western discourse shifting from governance/safety toward labor and ethics concerns. Occupation-based Reddit analysis showed IT professionals most positive (65.0%), Law most negative (44.0%), and Humanities/Social Sciences most neutral. Topic modeling revealed field-specific concerns — Artists on copyright/IP, Engineers on automation vs. productivity, Teachers on burnout/plagiarism, Law on e-discovery vs. legal judgment. 
🧑‍🏫 [Click to View Slides](assets/projects/Ewha2026Spring/DataEngineering/Prj1Slides.pdf)

### Regional Cultural Accessibility Clustering and Type-Specific Policy Design
*Term project for Management Science Theory, Ewha Womans University (Spring 2026)* 
South Korea's cultural welfare spending has steadily expanded physical infrastructure, yet access gaps persist even where facilities already exist: rural and low-income households, older adults, and people with disabilities still report far lower cultural participation than the national average. We argued this is because past policy allocated resources based on a single, undifferentiated supply metric (facility counts), which cannot distinguish a region that lacks infrastructure altogether from one that has facilities but lacks physical, informational, or economic access to them. Treating regional cultural accessibility as a multidimensional segmentation problem rather than a single ranked index, we asked whether unsupervised clustering could reveal distinct regional deprivation patterns that a composite score would blur together, and use that typology to design differentiated public-sector interventions. 
Using all 250 Korean municipalities (si-gun-gu) as the unit of analysis, we built a region-level dataset by integrating six public data sources — population statistics (KOSIS), barrier-free cultural/tourism sites (Korea Culture Information Service), Culture Card (문화누리카드) offline merchant listings (Korea Arts Council), performance venue records (KOPIS), sports facility data (Korea Sports Promotion Foundation), and cultural-facility usage indicators — into 42 population-normalized variables covering facility supply, barrier-free accessibility, voucher-merchant availability, and visitor activity. After comparing K-means and hierarchical clustering, we selected K-means with k = 5 (validated via the elbow method, a silhouette score of 0.53, and PCA visualization), yielding five interpretable regional types: infrastructure-rich hubs (25 regions), rural/underserved areas (75 regions), a single extreme outlier of facility hyper-concentration (Seoul's Jongno-gu), ordinary urban living-zone areas (141 regions, 56% of all municipalities), and tourism-driven but locally underserved areas (8 regions). For each type we diagnosed its dominant accessibility bottleneck and proposed a matching intervention — e.g., mobile cultural-service circuits for rural clusters versus voucher-merchant expansion for tourism-heavy but locally underserved areas — suggesting that data-driven regional typologies can support more targeted, equity-oriented cultural welfare budgeting than uniform, facility-count-based allocation. 
📄 [Click to View PDF (Korean Original)](assets/projects/Ewha2026Spring/BizDS/FinalReport)

### Finding Reliable Confidence Interval for the Odds Ratio through Simulation and Visualization
*Course project for Computational Statistics, Western Washington University (Winter 2025)*  
There are different ways to compute the confidence interval for the odds ratio (OR) with
common methods suggested by Woolf, Agresti, and Gart. This work focuses on finding the best
method for calculating the OR using simulation and visualization in different ways, focusing on
the analysis of the three common methods with and without Welch’s adjustment. Parameters
were chosen ranging from moderate values to extreme ones, with 1760 unique combinations in
total. Monte Carlo simulation was used for examination. According to the output, tw, ta, and zw
performed well with ta showing the most success. The only case in which ta failed was when
both population sample sizes were extremely small. In addition, further exploration was applied
to a real-world example of smoking status and lung cancer diagnosis. Other simulation studies
with 90% and 99% confidence levels were also performed.  
📄 [Click to View PDF](assets/projects/OddsRatio/MATH445_FinalReport.pdf)

### Relationship between Academic Achievement and Future Career and Living
*Course project for Calling BS in the Age of Big Data and AI, Western Washington University (Winter 2025)*  
We wanted to determine if certain education factors can help predict future careers and living outlooks through data analysis. Using the dataset 'Education & Career Success' sourced from Kaggle, we focused on two research questions, 'Which educational factors among university ranking, university GPA, and field of study are related to work-life balance?' and 'Is university GPA related to job offers, starting salary, and
years to promotion?' We used regression analysis, ANOVA, and decision tree for the analysis but we couldn't find a significant relationship between academic performance, future living, and success.  
🧑‍🏫 [Click to View Slides](assets/projects/AchievementFuture/CSCI497C_Final_Slides.pdf)  
📄 [Click to View PDF](assets/projects/AchievementFuture/CSCI497C_Final_Report.pdf)

### Using Three Tree-Based Ensemble Models to Assess the Risk of Myocardial Infarction
*Course project for the Data Structures, Ewha Womans University (Spring 2024)*  
Myocardial Infarction is caused by clotted blood block vessels. It brings sudden chest pain without early symptoms and has 33% death rate upon arrival at the hospital. So we built model for checking the risk of Myocardial Infarction easily and frequently. We used tree-based models because it is easy to understand the decision making process and the searching is fast in trees in terms of time complexity. We used Korea National Health and Nutrition Examination Survey (KNHANES) data for training the model. We trained decision tree, random forest, gradient boosting, and XGBoosting models and XGBoosting showed the best model performance with F1-score 0.9751.  
🧑‍🏫 [Click to View Slides](assets/projects/MyocardialInfarction/DataStructure_FinalProjectSlide.pdf)

### Building Drug Usage Prediction Model
*Course project for the Data Science Foundations, Ewha Womans University (Spring 2024)*  
Drug abuse is a huge social problem in the U.S. We tried to build prediction model to drug use frequency. We used National Survey on Drug Use and Health (NSDUH) data provided by Substance Abuse and Mental Health Service Administration (SAMHSA). The target variable we used was 'the number of days one used marijuana in one year.' We first trained a random forest regressor, whose cross-validated $R^2$ ranged from -6.8 to -2.7  and applied Yeo Johnson and Bayesian optimization for better performance which showed cross validated mean $R^2$ is 0.23. Also we further discussed about ethical issues involved in using data such as race, economic status, and educational level to predict drug use frequency.   
🧑‍🏫 [Click to View Slides (Korean Original)](assets/projects/Drug/DSFoundation_TermProject_slides.pdf)

### Logistic Regression Analysis of Lung Cancer Incidence by Smoking Duration and Occupation
*Course project for the Health Medical Database, Ewha Womans University (Spring 2024)*  
Lung cancer is one of the most common cancer in South Korea, and among men, it has the highest incidence rate. In addition, lung cancer has a poor prognosis. So I conducted logistic regression analysis of lung cancer incidence by smoking duration and occupation. I used Korea Health Panel Survey data for the analysis. The analysis showed that having an occupation with high asbestos exposure did not adequately explain the presence of lung cancer. In contrast, past smoking duration was a strong explanatory variable. We found that for every one-year increase in past smoking duration, the probability of not having lung cancer decreased by a factor of 0.9278, indicating that a longer history of smoking increases the likelihood of developing lung cancer. Regarding model fit, the R statistic suggested that the model was not adequate, whereas the Hosmer–Lemeshow test and the ROC curve’s AUC value indicated that the model had good fit.  
🧑‍🏫 [Click to View Slides (Korean Original)](assets/projects/LungCancer/BioMedicalDB_presentationSlides.pdf)  
📄 [Click to View PDF (Korean Original)](assets/projects/LungCancer/BioMedicalDB_FinalReport.pdf)









