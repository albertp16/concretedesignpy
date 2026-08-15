/* MIT License — Copyright (c) Albert Pamonag Engineering Consultancy
 *
 * Calculation-sheet renderer.
 * ---------------------------
 * Turns a report payload from /api/beam/flexure-report,
 * /api/beam/shear-report or /api/joint/shear-report into an A4 document.
 *
 * THIS FILE HOLDS NO ENGINEERING. It compares nothing, derives nothing and
 * rounds nothing that matters — every number, every verdict and every QAQC
 * comparison arrives already decided from the server. That separation is
 * the point: a renderer that recomputes is a second implementation, and
 * this repository has already shipped one of those by accident.
 *
 * UMD so the browser can script-tag it and Jest can require it.
 */
(function (root, factory) {
    if (typeof module === 'object' && module.exports) {
        module.exports = factory();
    } else {
        root.CalcSheet = factory();
    }
}(typeof self !== 'undefined' ? self : this, function () {
    'use strict';

    /* ── primitives ───────────────────────────────────────────────── */

    // Text from the server is trusted-but-escaped; text that originated in a
    // user's form field is not trusted at all. Escaping everything is
    // cheaper than remembering which is which.
    function esc(s) {
        if (s === null || s === undefined) { return ''; }
        return String(s)
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#39;');
    }

    // Server strings carry deliberate entities (&phi;, &sup2;, &mdash;) and
    // occasionally a <span class="dim">. Those are ours, so they pass; the
    // set of allowed tags is closed and tiny on purpose.
    //
    // The close tag is balanced against the open tag rather than allowed
    // outright. Allowing a bare </span> meant a REJECTED opening tag still
    // left its closing tag behind: '<span class="evil">x</span>' came out as
    // 'x</span>'. No attribute survives either way so it was never an
    // injection, but stray markup on a printed calculation sheet is a defect
    // in its own right, and the depth counter is two lines.
    var OPEN_OK = /^<span(?: class="(?:dim|ok|fail)")?>$/;
    var CLOSE = /^<\/span>$/;
    function trust(s) {
        if (s === null || s === undefined) { return ''; }
        var depth = 0;
        return String(s).replace(/<[^>]*>/g, function (tag) {
            if (OPEN_OK.test(tag)) { depth++; return tag; }
            if (CLOSE.test(tag) && depth > 0) { depth--; return tag; }
            return '';
        });
    }

    // Formatting is presentation only. `null` means the server declined to
    // produce a value — it is rendered as an em dash, never as 0.
    function num(v, dp) {
        if (v === null || v === undefined || v === '') { return '&mdash;'; }
        var n = Number(v);
        if (!isFinite(n)) { return '&mdash;'; }
        return n.toFixed(dp === undefined ? 2 : dp);
    }

    function attr(s) { return esc(s); }

    /* ── masthead ─────────────────────────────────────────────────── */

    function masthead(d, opts) {
        opts = opts || {};
        var meta = [];
        if (opts.project) { meta.push('<b>' + esc(opts.project) + '</b>'); }
        if (opts.member) { meta.push('Member: <b>' + esc(opts.member) + '</b>'); }
        if (opts.date) { meta.push(esc(opts.date)); }
        meta.push(esc((d.provenance && d.provenance.report_version) || ''));
        return '<div class="cs-masthead">'
             + '<div><h2>' + esc(d.title) + '</h2>'
             + '<p class="cs-basis">' + trust(d.basis) + '</p></div>'
             + '<div class="cs-meta">' + meta.join('<br>') + '</div>'
             + '</div>';
    }

    /* ── verdict ──────────────────────────────────────────────────── */

    // Three states, not two. `adequate: null` means no demand was supplied,
    // and a capacity with nothing to compare it against has no verdict —
    // rendering that as a fail would be as wrong as rendering it as a pass.
    function verdict(d) {
        var gov = d.governing_checks || [];
        var missing = (d.unavailable || []).length;
        var cls, head, body;

        if (d.adequate === null || d.adequate === undefined) {
            cls = 'none';
            head = 'Capacity reported &mdash; no demand supplied';
            body = 'No factored demand was given, so no demand/capacity '
                 + 'verdict is drawn. The capacities below stand on their own.';
        } else if (d.adequate) {
            cls = 'pass';
            head = 'All computed checks satisfied';
            body = 'Every check this sheet actually performed is satisfied.';
        } else {
            cls = 'fail';
            head = gov.length + ' check' + (gov.length === 1 ? '' : 's')
                 + ' not satisfied';
            body = gov.map(function (g) { return esc(g.check); }).join(' &middot; ');
        }

        var caveat = '';
        if (missing) {
            caveat = '<span class="cs-caveat"><b>' + missing + ' required check'
                   + (missing === 1 ? ' was' : 's were') + ' NOT performed by '
                   + 'this sheet</b> &mdash; see Checks Not Performed below. '
                   + '&ldquo;All computed checks satisfied&rdquo; is a '
                   + 'statement about the arithmetic, not about the member.'
                   + '</span>';
        }
        return '<div class="cs-verdict ' + cls + '"><strong>' + head
             + '</strong>' + body + caveat + '</div>';
    }

    /* ── the REFERENCES | CALCULATIONS | RESULT sheet ─────────────── */

    function step(s, alt) {
        var cls = alt ? ' cs-row-alt' : '';
        var status = s.status === 'ok' ? 'ok' : (s.status === 'fail' ? 'fail' : '');
        return '<div class="cs-ref' + cls + '">' + trust(s.ref) + '</div>'
             + '<div class="cs-calc' + cls + '">'
             + '<p class="cs-desc">' + trust(s.desc) + '</p>'
             + '<div class="cs-eq">' + trust(s.eq) + '</div></div>'
             + '<div class="cs-res' + cls + '">'
             + (status ? '<span class="' + status + '">' : '')
             + trust(s.result)
             + (status ? '</span>' : '')
             + '</div>';
    }

    function sheet(d) {
        var h = '<div class="cs-sheet">'
              + '<div class="cs-hd">References</div>'
              + '<div class="cs-hd">Calculations</div>'
              + '<div class="cs-hd">Result</div>';
        var alt = false;
        (d.sections || []).forEach(function (sec) {
            h += '<div class="cs-shead">' + esc(sec.heading) + '</div>';
            alt = false;
            (sec.steps || []).forEach(function (s) {
                h += step(s, alt);
                alt = !alt;
            });
        });
        return h + '</div>';
    }

    /* ── summary ──────────────────────────────────────────────────── */

    function summary(d) {
        var rows = (d.summary || []).map(function (it) {
            return '<div class="cs-item"><span class="cs-lab">'
                 + trust(it.label) + '</span><span class="cs-val">'
                 + trust(it.value)
                 + (it.note ? '<span class="cs-sub">' + trust(it.note) + '</span>' : '')
                 + '</span></div>';
        }).join('');
        return '<div class="cs-block"><h3>Results Summary</h3>'
             + '<div class="cs-summary">' + rows + '</div></div>';
    }

    /* ── QAQC ─────────────────────────────────────────────────────── */

    function qaqc(d) {
        var q = d.qaqc;
        if (!q || !q.checks || !q.checks.length) {
            return '<div class="cs-block"><h3>QAQC</h3>'
                 + '<p class="cs-empty">No QAQC block was returned. Do not '
                 + 'file this sheet.</p></div>';
        }
        var h = '<div class="cs-block"><h3>QAQC &mdash; Independent '
              + 'Recomputation (' + q.n_pass + ' of ' + q.n_total
              + ' pass)</h3>';
        h += '<p class="cs-note">' + esc(q.note) + ' Every expected and '
           + 'computed value below was derived on the server; this page '
           + 'compares nothing itself.</p>';
        if (!q.all_pass) {
            h += '<div class="cs-qaqc-fail">One or more QAQC checks FAILED. '
               + 'This report is not internally consistent &mdash; do not use '
               + 'these numbers, and report it as a software defect.</div>';
        }
        h += '<div class="cs-tablewrap">'
           + '<table class="cs-table"><thead><tr><th>Reported value</th>'
           + '<th>Independent method</th><th>Expected</th><th>Reported</th>'
           + '<th>Status</th></tr></thead><tbody>';
        q.checks.forEach(function (c) {
            h += '<tr><td>' + esc(c.name) + '</td>'
               + '<td>' + esc(c.method) + '</td>'
               + '<td class="num">' + num(c.expected, 4) + '</td>'
               + '<td class="num">' + num(c.computed, 4) + '</td>'
               + '<td class="mid"><span class="cs-chip ' + (c.pass ? 'pass' : 'fail')
               + '">' + (c.pass ? 'PASS' : 'FAIL') + '</span></td></tr>';
        });
        return h + '</tbody></table></div></div>';
    }

    /* ── advisories ───────────────────────────────────────────────── */

    // An advisory is part of the answer, so the full text is always in the
    // DOM. Nothing here truncates, filters or hides one.
    function advisories(d) {
        var list = d.advisories || [];
        if (!list.length) {
            return '<div class="cs-block"><h3>Advisories</h3>'
                 + '<p class="cs-empty">None returned.</p></div>';
        }
        var order = { critical: 0, warning: 1, info: 2 };
        var sorted = list.slice().sort(function (a, b) {
            return (order[a.severity] === undefined ? 3 : order[a.severity])
                 - (order[b.severity] === undefined ? 3 : order[b.severity]);
        });
        var h = '<div class="cs-block"><h3>Advisories (' + list.length
              + ')</h3><p class="cs-note">Engineering traps the arithmetic '
              + 'cannot rule out. These are part of the answer, not a '
              + 'footnote to it.</p>';
        sorted.forEach(function (a) {
            h += '<div class="cs-adv ' + attr(a.severity || 'info') + '">'
               + '<span class="cs-code">' + esc(a.code) + '</span>'
               + esc(a.text)
               + (a.clause ? '<span class="cs-clause">' + esc(a.clause) + '</span>' : '')
               + '</div>';
        });
        return h + '</div>';
    }

    /* ── checks not performed ─────────────────────────────────────── */

    function unavailable(d) {
        var list = d.unavailable || [];
        if (!list.length) {
            return '<div class="cs-block"><h3>Checks Not Performed</h3>'
                 + '<p class="cs-empty">None &mdash; every check this sheet '
                 + 'covers was performed.</p></div>';
        }
        var h = '<div class="cs-block"><h3>Checks Not Performed ('
              + list.length + ')</h3>'
              + '<p class="cs-note">This sheet has no model for the '
              + 'following. They are recorded rather than faked, and they '
              + 'remain the responsibility of the Engineer of Record.</p>'
              + '<div class="cs-tablewrap">'
              + '<table class="cs-table"><thead><tr><th>Check</th>'
              + '<th>Why</th><th>Clause</th></tr></thead><tbody>';
        list.forEach(function (u) {
            h += '<tr><td>' + esc(u.check) + '</td><td>' + esc(u.why)
               + '</td><td>' + esc(u.clause) + '</td></tr>';
        });
        return h + '</tbody></table></div></div>';
    }

    /* ── governing checks ─────────────────────────────────────────── */

    function governing(d) {
        var list = d.governing_checks || [];
        if (!list.length) { return ''; }
        var h = '<div class="cs-block"><h3>Governing Checks (' + list.length
              + ')</h3><div class="cs-tablewrap">'
              + '<table class="cs-table"><thead><tr><th>Check</th>'
              + '<th>D/C</th><th>Detail</th><th>Clause</th></tr></thead><tbody>';
        list.forEach(function (g) {
            h += '<tr><td>' + esc(g.check) + '</td>'
               + '<td class="num">' + num(g.dc, 3) + '</td>'
               + '<td>' + esc(g.detail) + '</td>'
               + '<td>' + esc(g.clause) + '</td></tr>';
        });
        return h + '</tbody></table></div></div>';
    }

    /* ── provenance ───────────────────────────────────────────────── */

    function provenance(d) {
        var p = d.provenance || {};
        return '<div class="cs-block"><h3>Provenance</h3>'
             + '<dl class="cs-prov">'
             + '<dt>Report version</dt><dd>' + esc(p.report_version) + '</dd>'
             + '<dt>Units</dt><dd>' + esc(p.units) + '</dd>'
             + '<dt>Code basis</dt><dd>' + esc(p.code_basis) + '</dd>'
             + '<dt>Clause register</dt><dd>' + esc(p.clause_register) + '</dd>'
             + '<dt>Disclaimer</dt><dd>' + esc(p.disclaimer) + '</dd>'
             + '</dl></div>';
    }

    /* ── the whole document ───────────────────────────────────────── */

    function render(d, opts) {
        if (!d) { return '<p class="cs-empty">No report to render.</p>'; }
        return '<div class="cs-doc">'
             + masthead(d, opts)
             + verdict(d)
             + summary(d)
             + sheet(d)
             + governing(d)
             + advisories(d)
             + unavailable(d)
             + qaqc(d)
             + provenance(d)
             + '</div>';
    }

    return {
        render: render,
        // exported for testing and reuse; each returns a fragment
        esc: esc,
        trust: trust,
        num: num,
        masthead: masthead,
        verdict: verdict,
        step: step,
        sheet: sheet,
        summary: summary,
        qaqc: qaqc,
        advisories: advisories,
        unavailable: unavailable,
        governing: governing,
        provenance: provenance
    };
}));
