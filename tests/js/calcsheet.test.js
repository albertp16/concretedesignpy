/* MIT License — Copyright (c) Albert Pamonag Engineering Consultancy
 *
 * Jest suite for the A4 calculation-sheet renderer.
 *
 * The renderer holds no engineering, so these tests hold none either. What
 * they assert is the set of properties that make a printed calculation sheet
 * trustworthy:
 *
 *   1. nothing the server sends is silently dropped -- an advisory, a
 *      not-performed check or a failed QAQC row must always reach the page
 *   2. a null is never rendered as a zero
 *   3. the three verdict states stay three
 *   4. untrusted text cannot inject markup
 */

const CalcSheet = require('../../concretedesignpy/webapp/static/js/calcsheet.js');
const { baseReport } = require('./fixtures');

function dom(html) {
    document.body.innerHTML = html;
    return document.body;
}

// ── primitives ────────────────────────────────────────────────────────

describe('esc', () => {
    test('escapes every character that can open a tag or attribute', () => {
        expect(CalcSheet.esc('<script>')).toBe('&lt;script&gt;');
        expect(CalcSheet.esc('a & b')).toBe('a &amp; b');
        expect(CalcSheet.esc('say "hi"')).toBe('say &quot;hi&quot;');
        expect(CalcSheet.esc("it's")).toBe('it&#39;s');
    });

    test('ampersand is escaped first, so entities are not double-decoded', () => {
        expect(CalcSheet.esc('&lt;')).toBe('&amp;lt;');
    });

    test('null and undefined render as empty, never as the words', () => {
        expect(CalcSheet.esc(null)).toBe('');
        expect(CalcSheet.esc(undefined)).toBe('');
    });

    test('zero survives -- it is a value, not an absence', () => {
        expect(CalcSheet.esc(0)).toBe('0');
    });
});

describe('trust', () => {
    test('keeps the closed set of span tags the server emits', () => {
        expect(CalcSheet.trust('<span class="dim">note</span>'))
            .toBe('<span class="dim">note</span>');
        expect(CalcSheet.trust('<span class="ok">SAFE</span>'))
            .toBe('<span class="ok">SAFE</span>');
    });

    test('strips anything outside that set', () => {
        expect(CalcSheet.trust('<script>alert(1)</script>')).toBe('alert(1)');
        expect(CalcSheet.trust('<img src=x onerror=1>')).toBe('');
        expect(CalcSheet.trust('<span class="evil">x</span>')).toBe('x');
        expect(CalcSheet.trust('<span onclick="x">y</span>')).toBe('y');
    });

    test('leaves the entities a calculation sheet is full of alone', () => {
        expect(CalcSheet.trust('&phi;Mn = 245 kN&middot;m'))
            .toBe('&phi;Mn = 245 kN&middot;m');
        expect(CalcSheet.trust('A_j mm&sup2;')).toBe('A_j mm&sup2;');
    });

    test('leaves MathJax delimiters intact', () => {
        const eq = '\\( \\beta_1 = 0.85 \\)';
        expect(CalcSheet.trust(eq)).toBe(eq);
    });
});

describe('num', () => {
    test('formats to the requested precision', () => {
        expect(CalcSheet.num(3.14159, 2)).toBe('3.14');
        expect(CalcSheet.num(3.14159, 4)).toBe('3.1416');
        expect(CalcSheet.num(245, 0)).toBe('245');
    });

    test('null is an em dash, NOT zero -- the whole point', () => {
        expect(CalcSheet.num(null, 2)).toBe('&mdash;');
        expect(CalcSheet.num(undefined, 2)).toBe('&mdash;');
        expect(CalcSheet.num('', 2)).toBe('&mdash;');
    });

    test('a real zero still prints as zero', () => {
        expect(CalcSheet.num(0, 2)).toBe('0.00');
    });

    test('non-finite values do not leak NaN or Infinity onto a sheet', () => {
        expect(CalcSheet.num(NaN, 2)).toBe('&mdash;');
        expect(CalcSheet.num(Infinity, 2)).toBe('&mdash;');
        expect(CalcSheet.num(-Infinity, 2)).toBe('&mdash;');
    });

    test('defaults to two decimal places', () => {
        expect(CalcSheet.num(1.239)).toBe('1.24');
    });
});

// ── masthead ──────────────────────────────────────────────────────────

describe('masthead', () => {
    test('carries the title, the basis and the report version', () => {
        const b = dom(CalcSheet.masthead(baseReport(), {}));
        expect(b.querySelector('h2').textContent).toBe('Beam Flexural Capacity');
        expect(b.querySelector('.cs-basis').textContent)
            .toContain('NSCP 2015 Section 422.2');
        expect(b.querySelector('.cs-meta').textContent)
            .toContain('beam-report-1.0');
    });

    test('project, member and date reach the sheet when supplied', () => {
        const b = dom(CalcSheet.masthead(baseReport(), {
            project: 'Palanca Building', member: 'B-12', date: '2026-08-15'
        }));
        const meta = b.querySelector('.cs-meta').textContent;
        expect(meta).toContain('Palanca Building');
        expect(meta).toContain('B-12');
        expect(meta).toContain('2026-08-15');
    });

    test('a project name cannot inject markup', () => {
        const b = dom(CalcSheet.masthead(baseReport(), {
            project: '<img src=x onerror=alert(1)>'
        }));
        expect(b.querySelector('img')).toBeNull();
        expect(b.querySelector('.cs-meta').textContent).toContain('<img');
    });
});

// ── verdict: three states ─────────────────────────────────────────────

describe('verdict', () => {
    test('adequate true reads as computed checks satisfied', () => {
        const b = dom(CalcSheet.verdict(baseReport({ adequate: true })));
        const v = b.querySelector('.cs-verdict');
        expect(v.classList.contains('pass')).toBe(true);
        expect(v.textContent).toContain('All computed checks satisfied');
    });

    test('adequate false names how many checks failed and which', () => {
        const b = dom(CalcSheet.verdict(baseReport({
            adequate: false,
            governing_checks: [
                { check: 'Flexural strength', dc: 1.4, detail: 'x', clause: 'y' },
                { check: 'Stirrup spacing', dc: 1.1, detail: 'x', clause: 'y' }
            ]
        })));
        const v = b.querySelector('.cs-verdict');
        expect(v.classList.contains('fail')).toBe(true);
        expect(v.textContent).toContain('2 checks not satisfied');
        expect(v.textContent).toContain('Flexural strength');
        expect(v.textContent).toContain('Stirrup spacing');
    });

    test('one failure is singular, because a sheet that says "1 checks" looks unread', () => {
        const b = dom(CalcSheet.verdict(baseReport({
            adequate: false,
            governing_checks: [{ check: 'Shear strength', dc: 2, detail: '', clause: '' }]
        })));
        expect(b.querySelector('.cs-verdict').textContent)
            .toContain('1 check not satisfied');
    });

    test('adequate null is its own state, neither pass nor fail', () => {
        const b = dom(CalcSheet.verdict(baseReport({
            adequate: null, has_demand: false
        })));
        const v = b.querySelector('.cs-verdict');
        expect(v.classList.contains('none')).toBe(true);
        expect(v.classList.contains('pass')).toBe(false);
        expect(v.classList.contains('fail')).toBe(false);
        expect(v.textContent).toContain('no demand supplied');
    });

    test('a pass verdict still carries the not-performed caveat', () => {
        const b = dom(CalcSheet.verdict(baseReport({ adequate: true })));
        const t = b.querySelector('.cs-verdict').textContent;
        expect(t).toContain('1 required check was NOT performed');
        expect(t).toContain('not about the member');
    });

    test('with nothing unavailable the caveat is absent', () => {
        const b = dom(CalcSheet.verdict(baseReport({
            adequate: true, unavailable: []
        })));
        expect(b.querySelector('.cs-caveat')).toBeNull();
    });

    test('the caveat pluralises', () => {
        const b = dom(CalcSheet.verdict(baseReport({
            unavailable: [
                { check: 'a', why: '', clause: '' },
                { check: 'b', why: '', clause: '' }
            ]
        })));
        expect(b.querySelector('.cs-verdict').textContent)
            .toContain('2 required checks were NOT performed');
    });
});

// ── the sheet grid ────────────────────────────────────────────────────

describe('sheet', () => {
    test('emits the three column headers once', () => {
        const b = dom(CalcSheet.sheet(baseReport()));
        const hd = [...b.querySelectorAll('.cs-hd')].map(e => e.textContent);
        expect(hd).toEqual(['References', 'Calculations', 'Result']);
    });

    test('every section heading and every step is rendered', () => {
        const d = baseReport();
        const b = dom(CalcSheet.sheet(d));
        expect(b.querySelectorAll('.cs-shead')).toHaveLength(d.sections.length);
        const steps = d.sections.reduce((n, s) => n + s.steps.length, 0);
        expect(b.querySelectorAll('.cs-ref')).toHaveLength(steps);
        expect(b.querySelectorAll('.cs-calc')).toHaveLength(steps);
        expect(b.querySelectorAll('.cs-res')).toHaveLength(steps);
    });

    test('each step keeps its reference, description, equation and result', () => {
        const b = dom(CalcSheet.sheet(baseReport()));
        expect(b.querySelectorAll('.cs-ref')[1].textContent)
            .toContain('Table 22.2.2.4.3, printed 438-439');
        expect(b.querySelectorAll('.cs-desc')[1].textContent)
            .toContain('Whitney stress-block factor');
        expect(b.querySelectorAll('.cs-eq')[1].textContent)
            .toContain('\\beta_1');
    });

    test('an ok status paints the result green, a fail red', () => {
        const ok = dom(CalcSheet.step(
            { ref: 'r', desc: 'd', eq: 'e', result: 'D/C = 0.8', status: 'ok' }, false));
        expect(ok.querySelector('.cs-res .ok')).not.toBeNull();
        const bad = dom(CalcSheet.step(
            { ref: 'r', desc: 'd', eq: 'e', result: 'D/C = 1.4', status: 'fail' }, false));
        expect(bad.querySelector('.cs-res .fail')).not.toBeNull();
    });

    test('a step with no status is painted neither way', () => {
        const b = dom(CalcSheet.step(
            { ref: 'r', desc: 'd', eq: 'e', result: '20 MPa', status: null }, false));
        expect(b.querySelector('.cs-res .ok')).toBeNull();
        expect(b.querySelector('.cs-res .fail')).toBeNull();
    });

    test('zebra striping alternates and restarts at each section', () => {
        const b = dom(CalcSheet.sheet(baseReport()));
        const refs = [...b.querySelectorAll('.cs-ref')];
        expect(refs[0].classList.contains('cs-row-alt')).toBe(false);
        expect(refs[1].classList.contains('cs-row-alt')).toBe(true);
        expect(refs[2].classList.contains('cs-row-alt')).toBe(false);
    });

    test('sections with no steps do not crash the renderer', () => {
        const b = dom(CalcSheet.sheet({ sections: [{ heading: 'Empty' }] }));
        expect(b.querySelector('.cs-shead').textContent).toBe('Empty');
    });

    test('no sections at all still yields a well-formed sheet', () => {
        const b = dom(CalcSheet.sheet({}));
        expect(b.querySelectorAll('.cs-hd')).toHaveLength(3);
    });
});

// ── summary ───────────────────────────────────────────────────────────

describe('summary', () => {
    test('renders one item per summary row, with its note', () => {
        const b = dom(CalcSheet.summary(baseReport()));
        expect(b.querySelectorAll('.cs-item')).toHaveLength(2);
        expect(b.querySelector('.cs-sub').textContent).toContain('0.9');
    });

    test('an item with no note omits the sub-line rather than printing empty', () => {
        const b = dom(CalcSheet.summary({
            summary: [{ label: 'c', value: '177 mm', note: '' }]
        }));
        expect(b.querySelector('.cs-sub')).toBeNull();
    });

    test('entities in labels and values survive', () => {
        const b = dom(CalcSheet.summary({
            summary: [{ label: 'Design moment, &phi;Mn', value: '245 kN&middot;m', note: '' }]
        }));
        expect(b.querySelector('.cs-lab').textContent).toContain('φMn');
        expect(b.querySelector('.cs-val').textContent).toContain('·');
    });
});

// ── QAQC ──────────────────────────────────────────────────────────────

describe('qaqc', () => {
    test('renders one row per check with expected, reported and a chip', () => {
        const b = dom(CalcSheet.qaqc(baseReport()));
        expect(b.querySelectorAll('tbody tr')).toHaveLength(2);
        expect(b.querySelectorAll('.cs-chip.pass')).toHaveLength(2);
        expect(b.querySelector('h3').textContent).toContain('2 of 2 pass');
    });

    test('a failed check is chipped FAIL and raises the banner', () => {
        const d = baseReport();
        d.qaqc.checks[1].pass = false;
        d.qaqc.n_pass = 1;
        d.qaqc.all_pass = false;
        const b = dom(CalcSheet.qaqc(d));
        expect(b.querySelectorAll('.cs-chip.fail')).toHaveLength(1);
        expect(b.querySelector('.cs-qaqc-fail')).not.toBeNull();
        expect(b.querySelector('.cs-qaqc-fail').textContent)
            .toContain('do not use these numbers');
    });

    test('no banner when everything passes', () => {
        expect(dom(CalcSheet.qaqc(baseReport())).querySelector('.cs-qaqc-fail'))
            .toBeNull();
    });

    test('a missing QAQC block is called out, not skipped silently', () => {
        const b = dom(CalcSheet.qaqc({}));
        expect(b.textContent).toContain('No QAQC block was returned');
        expect(b.textContent).toContain('Do not file this sheet');
    });

    test('an empty check list is treated the same way', () => {
        const b = dom(CalcSheet.qaqc({ qaqc: { checks: [] } }));
        expect(b.textContent).toContain('No QAQC block was returned');
    });

    test('a null expected value prints an em dash, not zero', () => {
        const d = baseReport();
        d.qaqc.checks[0].expected = null;
        d.qaqc.checks[0].computed = null;
        const b = dom(CalcSheet.qaqc(d));
        const cells = [...b.querySelectorAll('tbody tr')[0].querySelectorAll('.num')];
        expect(cells[0].textContent).toBe('—');
        expect(cells[0].textContent).not.toBe('0.0000');
    });

    test('the page states that it compares nothing itself', () => {
        expect(dom(CalcSheet.qaqc(baseReport())).textContent)
            .toContain('this page compares nothing itself');
    });
});

// ── advisories ────────────────────────────────────────────────────────

describe('advisories', () => {
    test('every advisory is rendered -- none is truncated or filtered', () => {
        const d = baseReport();
        const b = dom(CalcSheet.advisories(d));
        expect(b.querySelectorAll('.cs-adv')).toHaveLength(d.advisories.length);
        d.advisories.forEach(a => {
            expect(b.textContent).toContain(a.text);
            expect(b.textContent).toContain(a.code);
            expect(b.textContent).toContain(a.clause);
        });
    });

    test('critical sorts above warning, warning above info', () => {
        const b = dom(CalcSheet.advisories(baseReport()));
        const sev = [...b.querySelectorAll('.cs-adv')]
            .map(e => e.className.replace('cs-adv ', ''));
        expect(sev).toEqual(['critical', 'warning', 'info']);
    });

    test('an unknown severity sorts last rather than throwing', () => {
        const b = dom(CalcSheet.advisories({
            advisories: [
                { code: 'X', severity: 'banana', text: 't', clause: 'c' },
                { code: 'Y', severity: 'critical', text: 't', clause: 'c' }
            ]
        }));
        const codes = [...b.querySelectorAll('.cs-code')].map(e => e.textContent);
        expect(codes).toEqual(['Y', 'X']);
    });

    test('the count in the heading matches the number rendered', () => {
        const b = dom(CalcSheet.advisories(baseReport()));
        expect(b.querySelector('h3').textContent).toContain('(3)');
        expect(b.querySelectorAll('.cs-adv')).toHaveLength(3);
    });

    test('advisory text cannot inject markup', () => {
        const b = dom(CalcSheet.advisories({
            advisories: [{
                code: 'X', severity: 'info',
                text: '<script>alert(1)</script>', clause: 'c'
            }]
        }));
        expect(b.querySelector('script')).toBeNull();
        expect(b.textContent).toContain('<script>');
    });

    test('an empty list says so instead of rendering nothing', () => {
        const b = dom(CalcSheet.advisories({ advisories: [] }));
        expect(b.querySelector('.cs-empty').textContent).toContain('None returned');
    });

    test('a missing clause omits the line rather than printing undefined', () => {
        const b = dom(CalcSheet.advisories({
            advisories: [{ code: 'X', severity: 'info', text: 't' }]
        }));
        expect(b.querySelector('.cs-clause')).toBeNull();
        expect(b.textContent).not.toContain('undefined');
    });
});

// ── checks not performed ──────────────────────────────────────────────

describe('unavailable', () => {
    test('every not-performed check reaches the sheet with its clause', () => {
        const b = dom(CalcSheet.unavailable(baseReport()));
        expect(b.querySelectorAll('tbody tr')).toHaveLength(1);
        expect(b.textContent).toContain('As,min');
        expect(b.textContent).toContain('Section 9.6.1.2');
        expect(b.textContent).toContain('Engineer of Record');
    });

    test('an empty list is stated positively, not omitted', () => {
        const b = dom(CalcSheet.unavailable({ unavailable: [] }));
        expect(b.textContent).toContain('every check this sheet covers was performed');
    });
});

// ── governing checks ──────────────────────────────────────────────────

describe('governing', () => {
    test('renders nothing at all when there are none', () => {
        expect(CalcSheet.governing(baseReport())).toBe('');
    });

    test('renders each failing check with its D/C and clause', () => {
        const b = dom(CalcSheet.governing({
            governing_checks: [{
                check: 'Shear strength', dc: 1.42,
                detail: 'Vu exceeds phi Vn', clause: 'Section 22.5.1.1'
            }]
        }));
        expect(b.querySelector('h3').textContent).toContain('(1)');
        expect(b.textContent).toContain('Shear strength');
        expect(b.querySelector('.num').textContent).toBe('1.420');
        expect(b.textContent).toContain('Section 22.5.1.1');
    });

    test('a null D/C prints an em dash', () => {
        const b = dom(CalcSheet.governing({
            governing_checks: [{ check: 'x', dc: null, detail: '', clause: '' }]
        }));
        expect(b.querySelector('.num').textContent).toBe('—');
    });
});

// ── provenance ────────────────────────────────────────────────────────

describe('provenance', () => {
    test('carries version, units, basis, register and disclaimer', () => {
        const b = dom(CalcSheet.provenance(baseReport()));
        const t = b.textContent;
        expect(t).toContain('beam-report-1.0');
        expect(t).toContain('kN and kN.m');
        expect(t).toContain('NSCP 2015');
        expect(t).toContain('CLAUSES.md');
        expect(t).toContain('Implements code equations only');
    });

    test('a report with no provenance does not throw', () => {
        expect(() => dom(CalcSheet.provenance({}))).not.toThrow();
    });
});

// ── the whole document ────────────────────────────────────────────────

describe('render', () => {
    test('assembles every block in reading order', () => {
        const b = dom(CalcSheet.render(baseReport(), { project: 'P' }));
        const heads = [...b.querySelectorAll('.cs-block h3')].map(e => e.textContent);
        expect(heads[0]).toContain('Results Summary');
        expect(heads).toEqual(expect.arrayContaining([
            expect.stringContaining('Advisories'),
            expect.stringContaining('Checks Not Performed'),
            expect.stringContaining('QAQC'),
            expect.stringContaining('Provenance')
        ]));
    });

    test('the summary precedes the sheet, and QAQC follows both', () => {
        const html = CalcSheet.render(baseReport(), {});
        expect(html.indexOf('Results Summary')).toBeLessThan(html.indexOf('cs-sheet'));
        expect(html.indexOf('cs-sheet')).toBeLessThan(html.indexOf('QAQC'));
    });

    test('wraps everything in a single A4 document root', () => {
        const b = dom(CalcSheet.render(baseReport(), {}));
        expect(b.querySelectorAll('.cs-doc')).toHaveLength(1);
        expect(b.firstElementChild.className).toBe('cs-doc');
    });

    test('a null report says so instead of throwing', () => {
        expect(dom(CalcSheet.render(null)).textContent)
            .toContain('No report to render');
    });

    test('renders a shear-shaped report as readily as a flexure one', () => {
        const b = dom(CalcSheet.render(baseReport({
            title: 'Beam Shear Design',
            basis: 'NSCP 2015 Section 422.5 (= ACI 318M-14)',
            adequate: false,
            governing_checks: [{
                check: 'Section 22.5.1.2 cross-sectional limit',
                dc: 1.9, detail: 'Enlarge the section',
                clause: 'ACI 318-25M Section 22.5.1.2, printed 442'
            }]
        }), {}));
        expect(b.querySelector('h2').textContent).toBe('Beam Shear Design');
        expect(b.querySelector('.cs-verdict').classList.contains('fail')).toBe(true);
        expect(b.textContent).toContain('Enlarge the section');
    });

    test('produces no undefined or NaN anywhere on a complete sheet', () => {
        const t = dom(CalcSheet.render(baseReport(), { project: 'P' })).textContent;
        expect(t).not.toContain('undefined');
        expect(t).not.toContain('NaN');
        expect(t).not.toContain('[object Object]');
    });

    test('a sparse report still renders every block without throwing', () => {
        const b = dom(CalcSheet.render({ title: 'Bare', basis: '' }, {}));
        expect(b.querySelector('h2').textContent).toBe('Bare');
        expect(b.textContent).toContain('No QAQC block was returned');
        expect(b.textContent).toContain('None returned');
    });
});

// ── regression: span balancing ────────────────────────────────────────

describe('trust — span balancing', () => {
    test('a rejected opening tag takes its closing tag with it', () => {
        expect(CalcSheet.trust('<span class="evil">x</span>')).toBe('x');
        expect(CalcSheet.trust('<span onclick="x">y</span>')).toBe('y');
    });

    test('an unmatched closing tag is dropped', () => {
        expect(CalcSheet.trust('done</span>')).toBe('done');
        expect(CalcSheet.trust('</span></span>')).toBe('');
    });

    test('nested allowed spans still balance', () => {
        expect(CalcSheet.trust('<span class="ok">a<span class="dim">b</span></span>'))
            .toBe('<span class="ok">a<span class="dim">b</span></span>');
    });

    test('an allowed span keeps its closing tag when a rejected one precedes it', () => {
        expect(CalcSheet.trust('<div>x</div><span class="dim">y</span>'))
            .toBe('x<span class="dim">y</span>');
    });

    test('the rendered sheet never emits an unbalanced span', () => {
        const html = CalcSheet.render(baseReport(), { project: 'P' });
        const opens = (html.match(/<span\b/g) || []).length;
        const closes = (html.match(/<\/span>/g) || []).length;
        expect(opens).toBe(closes);
    });
});

// ── responsive: wide tables must not push the page sideways ───────────

describe('table overflow containment', () => {
    test('each data table is wrapped so it scrolls inside its own block', () => {
        const d = baseReport({
            governing_checks: [{ check: 'x', dc: 1.2, detail: 'y', clause: 'z' }]
        });
        const b = dom(CalcSheet.render(d, {}));
        const tables = b.querySelectorAll('.cs-table');
        expect(tables.length).toBe(3);          // qaqc, unavailable, governing
        tables.forEach(t => {
            expect(t.parentElement.className).toBe('cs-tablewrap');
        });
    });

    test('the rendered document has balanced divs', () => {
        const html = CalcSheet.render(baseReport(), { project: 'P' });
        const open = (html.match(/<div\b/g) || []).length;
        const close = (html.match(/<\/div>/g) || []).length;
        expect(open).toBe(close);
    });
});
