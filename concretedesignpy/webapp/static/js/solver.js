/**
 * ConcreteDesignPy — solver shell runtime.
 *
 * Powers the VIKTOR-style solver layout used by every calculator page:
 *
 *   - builds the view tab bar from `[data-view]` panels inside `.solver-views`
 *   - resizes Plotly charts when a hidden panel becomes visible
 *   - tracks run state (idle / stale / running / ok / error) by instrumenting
 *     window.fetch for POSTs to /api/*
 *   - optional auto-update: debounced form re-submit on parameter change
 *   - Ctrl/Cmd+Enter runs the solver from anywhere on the page
 *   - collapsible parametrization sections, mobile sidebar toggle
 *
 * Pure presentation glue: nothing in here computes, interprets or rewrites
 * results — pages keep their own payload builders and renderers.
 */
(function () {
    'use strict';

    var form = null;           // #calc-form, if the page has one
    var statusEl = null;       // .run-status chip
    var emptyEl = null;        // .solver-empty placeholder
    var autoTimer = null;
    var pendingRuns = 0;
    var runStart = 0;

    // ── Run status chip ──────────────────────────────────────────────

    function setStatus(state, text) {
        if (!statusEl) { return; }
        statusEl.dataset.state = state;
        statusEl.textContent = text;
    }

    // Instrument fetch so every page gets status + timing for free.
    var origFetch = window.fetch;
    window.fetch = function (url, opts) {
        var isApi = typeof url === 'string' && url.indexOf('/api/') === 0
                    && opts && (opts.method || 'GET').toUpperCase() === 'POST';
        if (isApi) {
            if (pendingRuns === 0) { runStart = performance.now(); }
            pendingRuns += 1;
            setStatus('running', 'Running…');
        }
        var p = origFetch.apply(this, arguments);
        if (!isApi) { return p; }

        return p.then(function (res) {
            // Peek at the JSON envelope without consuming the caller's body.
            res.clone().json().then(function (data) {
                finishRun(data && data.status === 'error'
                    ? (data.message || 'Calculation error') : null);
            }).catch(function () { finishRun(res.ok ? null : 'HTTP ' + res.status); });
            return res;
        }).catch(function (err) {
            finishRun(err && err.message ? err.message : 'Request failed');
            throw err;
        });
    };

    function finishRun(errMessage) {
        pendingRuns = Math.max(0, pendingRuns - 1);
        if (pendingRuns > 0) { return; }
        if (errMessage) {
            setStatus('error', 'Error — see panel');
        } else {
            var secs = ((performance.now() - runStart) / 1000).toFixed(1);
            setStatus('ok', 'Up to date · ' + secs + 's');
            hideEmptyState();
        }
    }

    function hideEmptyState() {
        if (emptyEl) { emptyEl.style.display = 'none'; }
    }

    // ── View tabs ────────────────────────────────────────────────────

    function resizePlots(panel) {
        if (!window.Plotly) { return; }
        var plots = panel.querySelectorAll('.js-plotly-plot');
        for (var i = 0; i < plots.length; i++) {
            try { window.Plotly.Plots.resize(plots[i]); } catch (e) { /* not mounted yet */ }
        }
    }

    function initTabs() {
        var views = document.querySelector('.solver-views');
        if (!views) { return; }
        var panels = views.querySelectorAll('[data-view]');
        if (!panels.length) { return; }

        var bar = views.querySelector('.view-tabs');
        if (!bar) {
            bar = document.createElement('nav');
            bar.className = 'view-tabs';
            views.insertBefore(bar, views.firstChild);
        }

        var tabs = [];
        Array.prototype.forEach.call(panels, function (panel, idx) {
            var tab = document.createElement('button');
            tab.type = 'button';
            tab.className = 'view-tab';
            tab.textContent = panel.getAttribute('data-view');
            tab.addEventListener('click', function () { activate(idx); });
            bar.appendChild(tab);
            tabs.push(tab);
        });

        function activate(idx) {
            Array.prototype.forEach.call(panels, function (panel, i) {
                var on = (i === idx);
                panel.hidden = !on;
                tabs[i].classList.toggle('active', on);
                if (on) {
                    resizePlots(panel);
                    if (window.MathJax && MathJax.typesetPromise) {
                        MathJax.typesetPromise([panel]).catch(function () {});
                    }
                }
            });
        }
        activate(0);

        // Pages may ask for a tab programmatically, e.g. after building the
        // calculation sheet: window.SolverShell.showView('Calc Sheet')
        window.SolverShell = window.SolverShell || {};
        window.SolverShell.showView = function (name) {
            Array.prototype.forEach.call(panels, function (panel, i) {
                if (panel.getAttribute('data-view') === name) { activate(i); }
            });
        };
    }

    // ── Auto-update + dirty tracking ─────────────────────────────────

    function autoKey() { return 'cdp-auto:' + location.pathname; }

    function initRunbar() {
        form = document.getElementById('calc-form');
        statusEl = document.querySelector('.run-status');
        emptyEl = document.querySelector('.solver-empty');
        if (!form) { return; }

        var autoBox = document.getElementById('solver-auto');
        if (autoBox) {
            autoBox.checked = localStorage.getItem(autoKey()) === '1';
            autoBox.addEventListener('change', function () {
                localStorage.setItem(autoKey(), autoBox.checked ? '1' : '0');
                if (autoBox.checked) { scheduleAutoRun(); }
            });
        }

        function onParamChange() {
            if (statusEl && statusEl.dataset.state !== 'running') {
                setStatus('stale', 'Parameters changed');
            }
            if (autoBox && autoBox.checked) { scheduleAutoRun(); }
        }

        form.addEventListener('input', onParamChange);
        form.addEventListener('change', onParamChange);

        function scheduleAutoRun() {
            clearTimeout(autoTimer);
            autoTimer = setTimeout(function () {
                if (typeof form.requestSubmit === 'function') { form.requestSubmit(); }
                else { form.dispatchEvent(new Event('submit', { cancelable: true })); }
            }, 700);
        }

        // Ctrl/Cmd+Enter runs from anywhere (including inside Handsontable).
        document.addEventListener('keydown', function (e) {
            if ((e.metaKey || e.ctrlKey) && e.key === 'Enter') {
                e.preventDefault();
                if (typeof form.requestSubmit === 'function') { form.requestSubmit(); }
            }
        });
    }

    // ── Collapsible parametrization sections ─────────────────────────

    function initCollapsibles() {
        var heads = document.querySelectorAll('.solver-params .form-section > h3');
        Array.prototype.forEach.call(heads, function (h3) {
            h3.addEventListener('click', function (e) {
                // Clicks on controls living inside the heading (e.g. the
                // rebar mode radio toggle) must not collapse the section.
                if (e.target !== h3 && e.target.closest('label,input,select,button')) { return; }
                h3.parentElement.classList.toggle('collapsed');
            });
        });
    }

    // ── Mobile sidebar ───────────────────────────────────────────────

    function initSidebar() {
        var toggle = document.getElementById('sidebar-toggle');
        var backdrop = document.getElementById('sidebar-backdrop');
        if (toggle) {
            toggle.addEventListener('click', function () {
                document.body.classList.toggle('sidebar-open');
            });
        }
        if (backdrop) {
            backdrop.addEventListener('click', function () {
                document.body.classList.remove('sidebar-open');
            });
        }
    }

    // ── Boot ─────────────────────────────────────────────────────────

    function boot() {
        initSidebar();
        initTabs();
        initRunbar();
        initCollapsibles();
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', boot);
    } else {
        boot();
    }
}());
