""" Builds the current version of the treatment parameters table (v3) for the recalibration; earlier versions, with the user's comments, are in versions/ """
import sys
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

out = sys.argv[1]
fig_qf = sys.argv[2]
version = 'v3'
na = '—'
body = 9.5 # Usable page width, inches

# Parameters table
cols = ['Symbol', 'Name', 'Description', 'Period', 'Current value', 'Source of current value', 'Suggested sources', 'Notes', 'Your input']
widths = [0.6, 1.0, 1.35, 0.85, 1.0, 1.0, 1.45, 1.75, 0.95]
widths = [w*body/sum(widths) for w in widths]

sections = []

sections.append(('A. Open: inputs and priors to set', 'New, or still open.', [
    ['α', 'Adherence factor', 'Scales the cure share for people who complete the course: cure-track share = α × c(t). Lets a regimen that is easier to adhere to raise cure.',
     '1995 on (one value)', na, na,
     'Imperial et al. 2018 (Nat Med): outcomes by doses missed, in pooled phase 3 trials.',
     ['Calibrated, with WHO\'s treatment success as a target (section B). Suggested prior U(0.7, 1.0).',
      'Success excludes relapses, so it separates α from p_relapse. But it pins only (1 − κ) × (1 − α × c), so α is identified only as well as κ is known.']],
    ['k_tiw', 'Thrice-weekly efficacy factor', 'Cure share of the thrice-weekly regimen relative to daily: c(t) = k_tiw × 0.95 before 2017.',
     '1995–2016, ramping to 1 over 2017', na, na,
     ['✓ Public-sector regimens (web search): thrice-weekly under RNTCP from 1997 (whole country by March 2006), as separate tablets in patient-wise boxes; daily dosing and fixed-dose combinations together from 2017 (ordered January 2017, national by about November 2017).',
      '✓ Indian studies cited at the 2017 rollout: relapse about 10% or more on intermittent vs about 5% on daily.',
      'Menzies et al. 2009 (PLoS Med): intermittent vs daily rifampicin.'],
     ['Suggested prior U(0.85, 1.0).',
      'The evidence is mostly on relapse; putting it on c makes the extra non-cures early failures (back at once) rather than relapses (back after a delay).',
      'WHO\'s recorded failure about halves across 2017 (1.0–1.5% of evaluated new cases in 2004–2015, 0.5–0.9% in 2017–2023). With κ constant, that fits non-cure among completers halving, i.e. k_tiw ≈ 0.95. Testing changed too (more molecular tests), so this is suggestive.',
      'Before 2006, districts outside RNTCP used the older National TB Programme regimens; treating them as thrice-weekly is a simplification. HIV-TB patients may have moved to daily dosing earlier (unverified).']],
    ['κ', 'Failures recorded as success', 'Share of people who complete a failing course that WHO\'s data record as treatment success ("treatment completed" without a test showing failure).',
     '1995 on (one value)', na, na,
     ['India TB Reports: share of new cases bacteriologically confirmed; share with an end-of-treatment follow-up test.',
      'Studies of culture at the end of treatment among patients recorded as completed (Indian cohorts).'],
     ['Calibrated. Suggested prior U(0.6, 0.95), or built as κ = 1 − b × t × s: b the bacteriologically confirmed share, t the share tested at the end, s the test\'s sensitivity for failure.',
      'The data need κ high: at α = 1, c = 0.95, recorded failure of 0.5–1.4% implies κ ≈ 0.7–0.9.',
      'A tight prior on κ is what lets the success target identify α.']],
    ['q(f)', 'Partial benefit after loss to follow-up', 'Chance that someone on the cure track lost after completing a fraction f of the course is cured, relative to completing it.',
     'All years', 'q(f) = f', 'Exp 35 placeholder',
     ['Truncated-regimen trials: BMRC 2-, 3- and 4-month regimens (Fox, Ellard & Mitchison 1999).',
      'Failed 4-month trials (REMoxTB, OFLOTUB, RIFAQUIN): relapse at f ≈ 2/3.',
      'Chang et al. (Hong Kong): relapse by doses taken.'],
     ['Suggested q(f) = f^a with prior a ~ U(0.3, 1). Figure below.',
      'Non-cures after partial treatment: suggest all failures for now (decision 1).']],
    ['ρ (or ω)', 'Early deaths the regimen can\'t prevent', 'TB death hazard on the cure track in the first phase (2 months), as a multiple of sym_dead; on the fail track such deaths are within μ_f. (ω is the share of deaths on treatment on the cure track.)',
     'All years', 'ω = 0 and 0.5 (exp 35\'s arms)', 'Exp 35 assumption',
     ['Waitt & Squire 2011 (IJTLD): deaths during and after treatment.',
      'Indian cohort studies of when deaths on treatment happen.'],
     ['Suggest ρ, with its prior set so that ω is about 0.25–0.75 at 2015 values: with ω fixed, a regimen that raises cure would also cut deaths on the cure track, which it shouldn\'t (decision 2).']],
]))

sections.append(('B. Settled: from data', 'WHO\'s treatment outcomes for India, all new cases, as shares of evaluated cases (exp 35\'s data).', [
    ['d(t)', 'Died share', 'Share of notified starts who die on treatment, from any cause (the model\'s TB and background deaths together).',
     '1995 on', '3.9%, constant', 'WHO outcome data; mean of 1999–2023 (exp 35)', na,
     ['Applies to every course in the model (no unnotified sector in this calibration). Equation below.',
      'The year-to-year changes (SD 0.55 points) look like noise. 1999–2003 from smear-positive cohorts only, scaled by × 0.85.']],
    ['ℓ(t)', 'Lost share', 'Share of starts lost to follow-up.',
     '1995–1999', '20.6%', 'WHO outcome data; least-squares knot at 1999 (exp 35)', na,
     ['Held at the 1999 value before 1999 (under 30% of each cohort evaluated during the rollout of DOTS). 1999–2003 from smear-positive cohorts only, scaled by × 1.05.',
      'For other countries, the knots are an input: choose them by piecewise-linear regression (e.g. the number of knots by BIC), in a general function.']],
    ['ℓ(t)', 'Lost share', 'As above.', '2003', '6.9%', 'As above', na, 'Linear between knots.'],
    ['ℓ(t)', 'Lost share', 'As above.', '2014', '5.6%', 'As above', na, 'Linear between knots.'],
    ['ℓ(t)', 'Lost share', 'As above.', '2023 on', '1.9%', 'As above', na, 'Constant after 2023. Residual SD of the fit 0.4 points.'],
    ['Success (target)', 'Treatment success', 'Cured plus completed, as a share of evaluated starts. A calibration target: matched to the model\'s success including κ (equation below).',
     '1999–2023', 'For example 0.88 (2005), 0.89 (2014), 0.94 (2023)', 'WHO outcome data: 1 − died − lost − failed', na,
     ['With d and ℓ matched exactly, the success target is the same as a target on recorded failure (1.0–1.5% in 2004–2015, 2.1% in 2016, 0.5–0.9% in 2017–2023).',
      'Which years and how to weight them: decision 5.']],
]))

sections.append(('C. Settled: derived, no input needed', 'Solved in each draw so that a course gives the died and lost shares exactly (equations below).', [
    ['h(t)', 'Loss-to-follow-up hazard', 'Rate of loss to follow-up, both tracks and every phase.',
     '1995–1999 / 2003 / 2014 / 2023 on', '0.47 / 0.15 / 0.12 / 0.04 per year', 'Solved from ℓ(t) (exp 35)', na,
     'Re-solved with α × c(t); values change little. Knots as for ℓ(t).'],
    ['μ_f(t)', 'Death hazard on the fail track', 'TB death rate on treatment for people whose regimen would fail; replaces sym_dead on treatment.',
     '1995 on', 'Median 0.20 per year (1995), 0.21 (2023)', 'Solved from d(t) (exp 35)', na, 'Moves with α, k_tiw and ρ in each draw.'],
    ['Background death', 'Deaths from other causes on treatment', 'WPP\'s rate by age and sex, as in every state.',
     'All years', 'WPP India', 'Demography (CONVENTIONS.md 10)', na, 'Counted within d(t).'],
]))

sections.append(('D. Settled: calibrated', 'Current values are exp 32\'s posterior: median (95% interval), IMIS seed 1 / seed 2. Linear between knots, constant after the last. All to be recalibrated.', [
    ['cov(t)', 'Annual treatment coverage', 'Share of people with symptomatic TB starting notified treatment within a year (rate −ln(1 − cov)).',
     'End of 1994', '0', 'Fixed (tb_LAI_TPT)', na,
     'To revisit later: treatment existed before 1995 (the National TB Programme from 1962), with low coverage and effectiveness.'],
    ['cov(t)', 'Annual treatment coverage', 'As above.', '2010',
     'Fraction of the 2024 value: 0.87 (0.54–0.99) / 0.80 (0.46–0.98); coverage ≈ 0.50 / 0.49', 'Calibrated; prior on the fraction U(0, 1)', na,
     ['Notified treatment only, since unnotified treatment is treated as having no benefit in this calibration (z = 1; section F).',
      'Other shapes over time: a later experiment.']],
    ['cov(t)', 'Annual treatment coverage', 'As above.', '2024', '0.57 (0.51–0.65) / 0.61 (0.53–0.71)', 'Calibrated; prior U(0.50, 0.95)', na, 'As above.'],
    ['p_relapse', 'Share of cures that relapse', 'Cured people who go to relapsing, then symptomatic TB.',
     '1995 on', '0.065 (0.030–0.092) / 0.073 (0.047–0.095)', 'Calibrated; prior U(0.02, 0.10)', na,
     'Prior U(0.02, 0.15): trials understate programme relapse (shorter follow-up, better adherence).'],
    ['over_detection', 'Over-detection of notifications', 'Clinically diagnosed notifications without TB: reported notifications = (1 + over_detection) × treatment starts.',
     '1995 on', '0.12 (0.05–0.18) / 0.14 (0.06–0.19)', 'Calibrated; prior U(0, 0.20)', na,
     'The upper bound is the share of notifications clinically diagnosed.'],
]))

sections.append(('E. Settled: fixed by design', 'CONVENTIONS.md 9.', [
    ['c(t)', 'Cure share of a completed regimen, full adherence', 'Chance that someone who completes the regimen, taking every dose, is cured.',
     '1995–2016', 'k_tiw × 0.95 (thrice-weekly)', 'Your input; k_tiw calibrated (section A)', na, 'See k_tiw.'],
    ['c(t)', 'Cure share of a completed regimen, full adherence', 'As above.', '2017', 'Linear ramp from k_tiw × 0.95 to 0.95', 'Rollout of daily dosing over 2017', na, na],
    ['c(t)', 'Cure share of a completed regimen, full adherence', 'As above.', '2018 on', '0.95 (daily)', 'Your input (fixed for now)',
     ['Control arms of phase 3 trials: REMoxTB (Gillespie et al. 2014), RIFAQUIN (Jindani et al. 2014), Study 31/A5349 (Dorman et al. 2021).'],
     'To be replaced by a narrow distribution from trials.'],
    ['D', 'Course length', 'Fixed course.', 'All years', '6 months', 'Regimen definition', na, 'Also for private courses, when they are added.'],
    ['w', 'Phase width', 'Course split into phases of fixed width (M = D/w), for hazards and partial benefit by phase.', 'All years', '2 months (3 phases)', 'Modelling choice (exp 35)', na, na],
    ['dur_relapse_tx', 'Mean time from cure to relapse', 'Exponential delay in relapsing.', 'All years', '1.5 years', 'tb_LAI_TPT', na, na],
    ['tx_from', 'States treatment starts from', 'Which states can start a course.', 'All years', 'Symptomatic TB only', 'Design (exp 35)', na, na],
]))

sections.append(('F. Later: unnotified private treatment', 'Not in this calibration (z = 1). The structure for when it is added; see "Private-sector treatment and notifications over time" below.', [
    ['z', 'Quality of unnotified treatment', 'Scales unnotified courses\' outcomes from those of notified courses (z = 0) to no benefit (z = 1): cure-track share (1 − z) × α × c(t), death hazard moved toward sym_dead by z. Same loss to follow-up.',
     'All years', '1 (no benefit)', 'Decided for this calibration',
     ['Udwadia et al. 2010 (PLoS One); Das et al. 2015 (Lancet Infect Dis): private prescribing quality.',
      'India TB Reports: outcomes of notified private patients.'],
     ['z = 1 is the same model as having no unnotified sector, except that people on a useless course would stop transmitting; leaving the sector out avoids that.',
      'Later: calibrated with a wide prior, or swept.']],
    ['u(t)', 'Unnotified share of courses', 'Share of all courses not notified: unnotified private treatment.',
     '1995–2014; 2024', na, na,
     ['Arinaminpathy et al. 2016 (Lancet Infect Dis): about 2.2 million private treatment episodes in 2014, from drug sales, against about 1.6 million notifications (to verify).',
      'Ni-kshay private notifications (India TB Reports).'],
     ['About 0.55 in 2014 if those numbers hold, held back to 1995; your suggestion 0 in 2024, linear between.',
      'Not needed while z = 1.']],
]))

# Rows the user has already accepted (in v1's comments, or in the conversation after v2), keyed by (symbol, period): shown in "Your input"
accepted = {
    ('d(t)', '1995 on'): 'Accepted in v1, with your change applied (share of notified starts).',
    ('ℓ(t)', '1995–1999'): 'Accepted in v1 (flexible knots for other countries noted).',
    ('ℓ(t)', '2003'): 'Accepted in v1.',
    ('ℓ(t)', '2014'): 'Accepted in v1.',
    ('ℓ(t)', '2023 on'): 'Accepted in v1.',
    ('h(t)', '1995–1999 / 2003 / 2014 / 2023 on'): 'Accepted in v1.',
    ('μ_f(t)', '1995 on'): 'Accepted in v1 (equation added below).',
    ('Background death', 'All years'): 'Accepted in v1.',
    ('cov(t)', 'End of 1994'): 'Accepted in v1 (flagged to revisit).',
    ('cov(t)', '2010'): 'Accepted in v1 (to recalibrate). Notified treatment only, from your z = 1 decision.',
    ('cov(t)', '2024'): 'Accepted in v1 (to recalibrate). Notified treatment only, from your z = 1 decision.',
    ('p_relapse', '1995 on'): 'Accepted in v1, with your prior U(0.02, 0.15).',
    ('over_detection', '1995 on'): 'Accepted in v1 (kept as is).',
    ('c(t)', '1995–2016'): 'Agreed after v2 (efficacy increase in 2017). Review k_tiw\'s prior in section A.',
    ('c(t)', '2017'): 'Agreed after v2.',
    ('c(t)', '2018 on'): 'Accepted in v1 (your value, 0.95).',
    ('D', 'All years'): 'Accepted in v1 (same for private courses).',
    ('w', 'All years'): 'Accepted in v1.',
    ('dur_relapse_tx', 'All years'): 'Accepted in v1.',
    ('tx_from', 'All years'): 'Accepted in v1.',
    ('z', 'All years'): 'Agreed after v2 (z = 1 for this calibration).',
    ('u(t)', '1995–2014; 2024'): 'Deferred with z = 1: no review needed now.',
}
for _, _, rows in sections:
    for row in rows:
        row.append(accepted.get((row[0], row[3]), ''))
assert sum(len(r) == 9 and r[8] != '' for _, _, rows in sections for r in rows) == len(accepted), 'An accepted row\'s key matches no row'

changes = [
    'Unnotified private treatment is left out of this calibration (z = 1, the quality of unnotified treatment, meaning no benefit), so coverage is notified treatment. Its structure for later (u(t) and z) is in section F; m_d is removed.',
    'The thrice-weekly regimen lowers c before 2017 through k_tiw (an efficacy factor), ramping to daily over 2017, instead of raising relapse. Regimen notes separate dosing frequency from formulation.',
    'WHO\'s treatment success is a calibration target, with κ, the share of failing completers recorded as success.',
    'Decided and removed from the decisions: approach 1 with z, z = 1 now; the thrice-weekly effect on efficacy; the unnotified share in 2024 (deferred with section F).',
]

decisions = [
    ['1', 'Non-cures after partial treatment', 'All failures (back at once), or some through relapsing (back after a delay).', 'All failures for now.'],
    ['2', 'Early deaths on treatment', 'ω (share of deaths on the cure track) or ρ (hazard relative to sym_dead).', 'ρ, with its prior set from a range of ω.'],
    ['3', 'IMIS weights', 'Fix exp 32\'s problem (badly fitting prior draws carrying large weights) in this calibration, or only report it.', 'Fix it.'],
    ['4', 'How progression curve B is counted (a target)', 'From transitions (ProgressionY1 0.060) or end-of-step snapshots, as tb_LAI_TPT (0.054).', 'Transitions, as for incidence.'],
    ['5', 'Success target: years and weight', 'Every year 1999–2023 with an overdispersed likelihood, or a few years (e.g. 2005, 2014, 2022).',
     'A few years, one on each side of 2017 at least, so 25 near-identical points don\'t outweigh the burden targets.'],
    ['6', 'Notifications over time as a target', 'Add with this calibration, or later with unnotified treatment (section F).', 'Later, with section F (below).'],
]

# Document
doc = Document()
sec = doc.sections[0]
sec.orientation = WD_ORIENT.LANDSCAPE
sec.page_width, sec.page_height = Inches(11), Inches(8.5)
for side in ['left_margin', 'right_margin', 'top_margin', 'bottom_margin']:
    setattr(sec, side, Inches(0.75))
doc.styles['Normal'].font.name = 'Calibri'
doc.styles['Normal'].font.size = Pt(10)
for name, size in [('Title', 18), ('Heading 1', 13), ('Heading 2', 11)]:
    s = doc.styles[name]
    s.font.name = 'Calibri'
    s.font.size = Pt(size)
    s.font.color.rgb = RGBColor(0x1F, 0x3A, 0x5F)

def shade(cell, fill):
    shd = OxmlElement('w:shd')
    shd.set(qn('w:val'), 'clear')
    shd.set(qn('w:color'), 'auto')
    shd.set(qn('w:fill'), fill)
    cell._tc.get_or_add_tcPr().append(shd)

def fill_cell(cell, value, size=8.5, bold=False):
    cell.text = ''
    for i, text in enumerate(value if isinstance(value, list) else [value]):
        p = cell.paragraphs[0] if i == 0 else cell.add_paragraph()
        p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.space_before = Pt(0)
        r = p.add_run(text)
        r.font.size = Pt(size)
        r.bold = bold

def make_table(header, rows_by_section, widths, input_col=None):
    t = doc.add_table(rows=1, cols=len(header))
    t.style = 'Table Grid'
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    for i, h in enumerate(header):
        fill_cell(t.rows[0].cells[i], h, bold=True)
        shade(t.rows[0].cells[i], 'D9E2F3')
    el = OxmlElement('w:tblHeader')
    el.set(qn('w:val'), 'true')
    t.rows[0]._tr.get_or_add_trPr().append(el)
    for title, note, rows in rows_by_section:
        if title:
            m = t.add_row().cells
            m = m[0].merge(m[-1])
            fill_cell(m, title + (f': {note}' if note else ''), size=9, bold=True)
            shade(m, 'EDEDED')
        for row in rows:
            r = t.add_row()
            for i, v in enumerate(row + [''] * (len(header) - len(row))):
                fill_cell(r.cells[i], v, bold=(i == 0))
            if input_col is not None:
                shade(r.cells[input_col], 'E2EFDA' if row[input_col:] and row[input_col] else 'FFF8DC') # Green: already accepted; yellow: to review
    layout = OxmlElement('w:tblLayout')
    layout.set(qn('w:type'), 'fixed')
    t._tbl.tblPr.append(layout)
    for i, gc in enumerate(t._tbl.tblGrid.findall(qn('w:gridCol'))):
        gc.set(qn('w:w'), str(int(widths[i]*1440)))
    for row in t.rows:
        k = 0
        for tc in row._tr.findall(qn('w:tc')): # Merged cells span several grid columns
            span = tc.tcPr.find(qn('w:gridSpan')) if tc.tcPr is not None else None
            n = int(span.get(qn('w:val'))) if span is not None else 1
            w = tc.get_or_add_tcPr().get_or_add_tcW()
            w.type, w.w = 'dxa', int(sum(widths[k:k + n])*1440)
            k += n
    return t

def para(text, bold_lead=None):
    p = doc.add_paragraph()
    if bold_lead:
        p.add_run(bold_lead).bold = True
    p.add_run(text)
    return p

def bullets(items):
    for item in items:
        p = doc.add_paragraph(style='List Bullet')
        if isinstance(item, tuple):
            p.add_run(item[0]).bold = True
            p.add_run(item[1])
        else:
            p.add_run(item)

def equation(text):
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Inches(0.4)
    r = p.add_run(text)
    r.font.name = 'Cambria Math'
    r.font.size = Pt(10.5)

doc.add_paragraph(f'Treatment parameters for the recalibration ({version})', style='Title')
para('Every treatment parameter in the planned recalibration with treatment outcomes: the inputs still open (section A), those settled from data, derived, calibrated or fixed by design (B–E), and unnotified private treatment, deferred (F). Please fill in the yellow cells of the "Your input" column; green cells are rows you have already accepted, which need no further review. Earlier versions, with your comments, are in versions/.')
para('Current values are from exp 35 (death and loss to follow-up on treatment) and exp 32\'s calibrated posterior. Sources marked ✓ were checked by web search; the others were written from my background knowledge and still need checking (authors, years and values).')
doc.add_heading('Changes from v2', level=2)
bullets(changes)

doc.add_heading('Parameters', level=1)
make_table(cols, sections, widths, input_col=8)

doc.add_heading('Open decisions', level=1)
dw = [0.3, 1.8, 3.4, 2.6, 1.4]
make_table(['#', 'Decision', 'Options', 'Suggestion', 'Your choice'], [(None, None, decisions)], [w*body/sum(dw) for w in dw], input_col=4)

doc.add_heading('Private-sector treatment and notifications over time', level=1)
para('The structure for unnotified treatment is approach 1: courses are notified or unnotified, an unnotified share u(t) shrinks over time, and unnotified courses have the same loss to follow-up but outcomes scaled by a quality z, from those of notified courses (z = 0) to no benefit (z = 1). This calibration fixes z = 1, which is the same model as exp 32\'s, with no unnotified sector.', bold_lead='Structure. ')
para('What z = 1 leaves out:', bold_lead='')
bullets([
    ('Natural history absorbs it. ', 'Whatever unnotified treatment really does to burden is attributed to TB death, self-cure and the duration of disease.'),
    ('Private-sector engagement. ', 'Notifying courses that were already happening would look like new treatment, overstating its impact.'),
    ('The rise in notifications. ', 'It is read as more treatment (below).'),
])
para('Notifications rise in the data to 173 per 100,000 in 2024, while the model\'s peak around 2009 and fall (exp 32). Much of the rise is private notification. With z = 1, every extra notified course is extra effective treatment, so fitting the rise would make burden fall faster: the same direction as the existing tension (mortality and prevalence already fall too fast). The rise can be fitted without speeding the decline only if unnotified treatment worked (z well below 1) and u fell, so that the extra notifications are mostly courses that were already happening. So notifications over time belong with section F, not with z = 1.', bold_lead='Fitting notifications over time. ')
para('Treatment starts could instead be forced from the data: starts = notifications ÷ (1 + over_detection), i.e. a per-person rate of starts ÷ the model\'s symptomatic TB each step. This replaces the coverage knots with data and matches notifications exactly. But:', bold_lead='Notifications as a direct input. ')
bullets([
    'It is a count per population, so the rate per person with symptomatic TB depends on the model\'s own prevalence. It runs out of people to treat where modelled prevalence is low, and needs a cap.',
    'With z = 1 it builds the rise in notifications in as new treatment, which forces the fast decline rather than testing it.',
    'Notifications must be defined as model starts are: the model counts every start, including restarts after failure, loss to follow-up and relapse, so the data series must include retreatment.',
    'The prevalence ÷ notifications target then constrains only prevalence.',
])
para('Keep coverage calibrated with the 2020 notifications target in this calibration. Add notifications over time with section F, and consider forcing then, for notified starts, with u and z calibrated.', bold_lead='Suggestion. ')

doc.add_heading('Outcome equations', level=1)
para('A course has M phases of width τ = D/M (here 3 of 2 months). It is on the cure track with probability π(t) = α × c(t), and on the fail track otherwise. In phase k on track j, TB death (μ_jk), background death (b, the age-weighted rate of people starting treatment) and loss to follow-up (h) compete, with total exit rate:')
equation('λ_jk = μ_jk + b + h')
para('The share of the track still on treatment at the start of phase k, and the shares dying and lost, are:')
equation('S_jk = exp(−τ × (λ_j1 + … + λ_j,k−1)),   S_j1 = 1')
equation('d = Σ_j π_j Σ_k S_jk × (μ_jk + b) / λ_jk × (1 − exp(−λ_jk τ))')
equation('ℓ = Σ_j π_j Σ_k S_jk × h / λ_jk × (1 − exp(−λ_jk τ))')
equation('cured = π × [ S_cure,M+1 + Σ_k S_cure,k × q_k × h / λ_cure,k × (1 − exp(−λ_cure,k τ)) ]')
para('Here π_cure = π and π_fail = 1 − π. The fail track has μ_fail,k = μ_f in every phase, which includes its early deaths. The cure track has the early-death hazard ρ × sym_dead in the first phase and no TB death after it. q_k = q(f_k) at the phase\'s midpoint f_k. Given d, ℓ, π, ρ and b, the equations for d and ℓ are solved for μ_f and h (tbode4.outcome_shares() and carve_out()).')
para('WHO\'s success counts completers, not people cured after loss to follow-up (WHO records them as lost), and includes failing completers not recorded as failures:')
equation('success = π × S_cure,M+1 + κ × (1 − π) × S_fail,M+1')
equation('recorded failure = (1 − κ) × (1 − π) × S_fail,M+1 = 1 − success − d − ℓ')
para('Since d and ℓ are matched exactly, the success target constrains (1 − κ) × (1 − π): α through π, given κ.')

doc.add_heading('Partial benefit after loss to follow-up: suggested shapes', level=1)
doc.add_picture(fig_qf, width=Inches(6.5))
para('q(f) = f^a for a across the suggested prior (0.3 to 1). The model applies q at the midpoint of each 2-month phase, so with a = 0.5 someone lost in months 0–2 is cured with probability 0.41 of a completer\'s, in months 2–4 0.71, and in months 4–6 0.91. (Script: make_qf_figure.py.)')

doc.save(out)
print('saved', out)
