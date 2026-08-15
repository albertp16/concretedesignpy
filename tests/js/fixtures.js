/* MIT License — Copyright (c) Albert Pamonag Engineering Consultancy
 *
 * Report fixtures for the Jest suite.
 *
 * These are trimmed but STRUCTURALLY REAL: every key the renderer touches is
 * present and has the shape the Python side actually emits. A fixture that
 * has drifted from the server is a test that passes while the page is
 * broken, so `tests/test_e2e_reports.py::test_fixture_shape_matches_the_server`
 * asserts the two agree.
 */

function baseReport(over) {
    var d = {
        provenance: {
            report_version: 'beam-report-1.0 (note 14 rectification, 2026-08-15)',
            units: 'mm, MPa internally; kN and kN.m in every reported value',
            code_basis: 'NSCP 2015 (= ACI 318M-14).',
            clause_register: 'CLAUSES.md',
            disclaimer: 'Implements code equations only.'
        },
        title: 'Beam Flexural Capacity',
        basis: 'NSCP 2015 Section 422.2 (= ACI 318M-14)',
        request_echo: { fc: 20, fy: 420, b: 250, h: 565 },
        sections: [
            {
                heading: '1. Design inputs',
                steps: [
                    {
                        ref: '--',
                        desc: 'Concrete cylinder strength',
                        eq: '\\( f\'_c = 20.0 \\text{ MPa} \\)',
                        result: '20.0 MPa',
                        status: null
                    },
                    {
                        ref: 'ACI 318-25M Table 22.2.2.4.3, printed 438-439',
                        desc: 'Whitney stress-block factor.',
                        eq: '\\( \\beta_1 = 0.85 \\)',
                        result: '0.85',
                        status: null
                    }
                ]
            },
            {
                heading: '2. Flexural capacity',
                steps: [
                    {
                        ref: 'ACI 318-25M Section 21.2.1, printed 429',
                        desc: 'Design flexural strength.',
                        eq: '\\( \\phi M_n = 245.25 \\)',
                        result: '&phi;Mn = 245.25 kN&middot;m',
                        status: null
                    },
                    {
                        ref: 'ACI 318-25M Section 9.5.1.1, printed 148',
                        desc: 'Flexural demand/capacity ratio.',
                        eq: '\\( M_u / \\phi M_n = 0.815 \\)',
                        result: 'D/C = 0.815',
                        status: 'ok'
                    }
                ]
            }
        ],
        summary: [
            { label: 'Neutral axis, c', value: '177.64 mm', note: '' },
            { label: 'Design moment, &phi;Mn', value: '245.25 kN&middot;m', note: '&phi; = 0.9' }
        ],
        qaqc: {
            checks: [
                {
                    name: 'beta1',
                    method: 'Table 22.2.2.4.3 longhand from f\'c',
                    expected: 0.85,
                    computed: 0.85,
                    pass: true
                },
                {
                    name: 'Design moment phi*Mn',
                    method: 'phi * Mn from the two reported values',
                    expected: 245.25,
                    computed: 245.25,
                    pass: true
                }
            ],
            n_pass: 2,
            n_total: 2,
            all_pass: true,
            note: 'Each row re-derives a reported value along an independent path.'
        },
        advisories: [
            {
                code: 'F-EDITION',
                severity: 'info',
                text: 'This sheet is NSCP 2015, not ACI 318-19.',
                clause: 'ACI 318-25M Table 21.2.2, printed 430'
            },
            {
                code: 'F-ASMIN',
                severity: 'critical',
                text: 'As,min is NOT checked by this sheet.',
                clause: 'ACI 318-25M Section 9.6.1.2, printed 149'
            },
            {
                code: 'F-RHO',
                severity: 'warning',
                text: 'The rho <= 0.025 cap is NOT checked.',
                clause: 'ACI 318-25M Section 18.6.3.1'
            }
        ],
        governing_checks: [],
        unavailable: [
            {
                check: 'As,min',
                why: 'not implemented',
                clause: 'ACI 318-25M Section 9.6.1.2, printed 149'
            }
        ],
        adequate: true,
        complete: false,
        has_demand: true
    };
    return Object.assign(d, over || {});
}

module.exports = { baseReport: baseReport };
