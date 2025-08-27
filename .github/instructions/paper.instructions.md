# Scientific Writing Guidelines (AWG-aligned)

These rules tell Copilot how to draft and revise scientific prose consistent with *Academic Writing for Graduate Students* (Swales & Feak) and common research-article norms (e.g., CARS introductions, IMRaD). When generating or editing text, follow ALL items below.

## 0) Integrity & boundaries
- **Never invent or fabricate references, data, or DOIs.** If a citation is missing, insert a placeholder like “[REF]” and ask for sources.  
- **Keep claims appropriately cautious (hedged) unless high-quality evidence warrants stronger phrasing.** Prefer “suggests/indicates/likely/may” over categorical claims when uncertainty exists. :contentReference[oaicite:1]{index=1}

## 1) Article macro-structure
- **Use IMRaD** for empirical papers (Introduction, Methods, Results, Discussion); add Abstract, Keywords, Acknowledgments, Funding, Data/Code availability as needed. :contentReference[oaicite:2]{index=2}
- **Craft Introductions via the CARS model** (three moves):  
  - **Move 1 – Establish territory:** define the research area and why it matters; synthesize key literature.  
  - **Move 2 – Establish a niche:** identify a gap/limitation, controversy, or unanswered question.  
  - **Move 3 – Occupy the niche:** state the study’s purpose, approach, and main contributions/hypotheses.  
  In Introductions, show stance with cautious language and properly integrated citations. :contentReference[oaicite:3]{index=3}

## 2) Paragraph & information flow
- **Old→New / Given→New:** start sentences/paragraphs with familiar information, end with new or heavier information; maintain theme–rheme progression for coherence. Use transitions (however, therefore, for example). :contentReference[oaicite:4]{index=4}
- **One purpose per paragraph:** topic sentence → development (evidence/reasoning) → wrap-up/transition.

## 3) Stance, hedging, and reporting others’ work
- **Hedge appropriately** (e.g., *may, might, appears, is likely, suggests, tends to*). Avoid over-certainty. Calibrate strength to evidence. :contentReference[oaicite:5]{index=5}
- **Use precise reporting verbs** to signal stance (e.g., *report, find, show, argue, contend, challenge*), and prefer present simple for general claims about prior literature; past for specific, completed studies. :contentReference[oaicite:6]{index=6}
- **Integrate sources** smoothly (not quote-dumping): summarize/paraphrase + citation; state how a cited finding relates to your research question.

## 4) Style & sentence-level choices
- **Favor strong verbs over nominalizations** (e.g., prefer *we analyzed* to *an analysis was conducted*), but allow nominalization when it improves cohesion or focus. Keep sentences clear and not overpacked. :contentReference[oaicite:7]{index=7}
- **Voice & tense defaults:**  
  - *Methods/Results* → past tense; mostly active unless agent is unimportant (passive acceptable).  
  - *General truths/lit review* → present.  
  - *Goals / paper organization* → present/future.
- **Terminology & definitions:** define key terms on first use; maintain consistent units, symbols, and abbreviations.

## 5) Section-specific expectations
- **Abstract (150–300 words unless journal dictates):** background (1–2 sentences) → objective → brief methods → **key numeric results** with effect sizes/uncertainty if available → conclusion/implication. Avoid citations and undefined acronyms.
- **Methods:** sufficient detail for replication; name protocols/software (with versions) and statistical approach.  
- **Results:** order sub-sections logically; **report numbers** (N, estimates, intervals, p-values as relevant); point to figures/tables.  
- **Discussion:** (i) recap key findings (no data dump), (ii) interpret in light of literature, (iii) limitations, (iv) implications and future work.
- **Figures/Tables:** each must be cited in text and have a **stand-alone** caption.

## 6) Rhetorical templates (use/adapt, don’t copy verbatim)
- **Move 1 (territory):** “Recent work has examined …; however, understanding of … remains limited.”  
- **Move 2 (niche):** “Nonetheless, prior studies have not addressed … / evidence is inconclusive regarding …”  
- **Move 3 (occupy):** “Here, we investigate … by … We show that …, which suggests …”
- **Hedged claims:** “These results **suggest/indicate** that …”; avoid “prove” unless logically warranted.

## 7) Ethical & transparency notes
- **State approvals/consent/permits** where applicable.  
- Include **data/code availability** statements and preregistration info where relevant.  
- **No hallucinated refs**: if unsure, write `[REF]` and request sources.

## 8) Revision checklist for Copilot (always run on drafts)
- Increase coherence with Given→New ordering.  
- Tighten nominalizations and prune redundancy.  
- Verify hedging is adequate; reduce booster words (e.g., *clearly, obviously*).  
- Check that Introduction follows CARS and Discussion includes limitations.  
- Ensure all claims are either supported by citations or flagged `[REF]`.

