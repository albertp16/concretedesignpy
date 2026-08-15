/**
 * ConcreteDesignPy Web App - Client-side JavaScript
 * Handles form submission and result rendering for all calculators.
 */

/**
 * Submit a calculation to the API and render results.
 * @param {string} url - API endpoint
 * @param {object} payload - JSON request body
 * @param {function} [chartCallback] - Optional callback for chart rendering
 */
function submitCalc(url, payload, chartCallback, errorCallback) {
    var resultsDiv = document.getElementById('results');
    var contentDiv = document.getElementById('results-content');

    if (resultsDiv) resultsDiv.style.display = 'block';
    if (contentDiv) contentDiv.innerHTML = '<div class="loading">Calculating...</div>';

    fetch(url, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
    })
    .then(function(res) { return res.json(); })
    .then(function(data) {
        if (data.status === 'error') {
            if (contentDiv) contentDiv.innerHTML = '<div class="error-msg">Error: ' + escapeHtml(data.message) + '</div>';
            if (errorCallback) errorCallback(data.message);
            return;
        }
        if (contentDiv) contentDiv.innerHTML = renderResult(data.result);
        if (chartCallback) {
            chartCallback(data.result);
        }
        // Run feedback lives in the solver status chip; a toast per run
        // would be noise. Toasts are reserved for errors.
    })
    .catch(function(err) {
        if (contentDiv) contentDiv.innerHTML = '<div class="error-msg">Request failed: ' + escapeHtml(err.message) + '</div>';
        if (errorCallback) errorCallback(err.message);
    });
}

/**
 * Render a result object as an HTML table.
 */
function renderResult(obj) {
    if (typeof obj !== 'object' || obj === null) {
        return '<p>' + escapeHtml(String(obj)) + '</p>';
    }

    var html = '<table class="result-table">';

    Object.keys(obj).forEach(function(key) {
        var val = obj[key];

        // Skip SVG/chart fields, large arrays, and detail objects from main table
        if (key === 'svg' || key === 'svg_pm' || key === 'svg_rebar' || key === 'svg_section' || key === 'chart_pm') return;
        if (key === 'rebar_layout' || key === 'rebar_forces') return;
        if (key === 'demand_check' || key === 'bar_depths' || key === 'max_moment_point') return;
        if (key === 'balanced_steel_detail' || key === 'bar_areas') return;
        if (key === 'geometry') return;
        if (Array.isArray(val) && val.length > 20) {
            return;
        }

        html += '<tr>';
        html += '<th>' + formatKey(key) + '</th>';

        if (Array.isArray(val)) {
            if (val.length > 0 && typeof val[0] === 'object') {
                html += '<td>' + renderArrayTable(val) + '</td>';
            } else {
                html += '<td>' + escapeHtml(JSON.stringify(val)) + '</td>';
            }
        } else if (typeof val === 'object' && val !== null) {
            html += '<td>' + renderResult(val) + '</td>';
        } else {
            html += '<td>' + formatValue(key, val) + '</td>';
        }

        html += '</tr>';
    });

    html += '</table>';
    return html;
}

/**
 * Render an array of objects as a sub-table.
 */
function renderArrayTable(arr) {
    if (arr.length === 0) return '(empty)';

    var keys = Object.keys(arr[0]);
    var html = '<table class="result-table" style="font-size:.8rem">';
    html += '<tr>';
    keys.forEach(function(k) { html += '<th>' + formatKey(k) + '</th>'; });
    html += '</tr>';

    arr.forEach(function(row) {
        html += '<tr>';
        keys.forEach(function(k) {
            html += '<td>' + formatValue(k, row[k]) + '</td>';
        });
        html += '</tr>';
    });

    html += '</table>';
    return html;
}

/**
 * Format a key name for display.
 */
function formatKey(key) {
    return escapeHtml(
        key.replace(/_/g, ' ')
           .replace(/\b\w/g, function(c) { return c.toUpperCase(); })
    );
}

/**
 * Format a value with status highlighting.
 */
function formatValue(key, val) {
    var str = String(val);
    var lower = str.toLowerCase();

    // Status indicators
    if (lower === 'pass' || lower === 'ok' || lower === 'compliant'
        || lower === 'tension-controlled' || lower === 'may neglect') {
        return '<span class="status-pass">' + escapeHtml(str) + '</span>';
    }
    if (lower === 'fail' || lower === 'noncompliant'
        || lower === 'reinforcement required' || lower === 'exceeds limit') {
        return '<span class="status-fail">' + escapeHtml(str) + '</span>';
    }
    if (lower === 'must consider' || lower === 'transition'
        || lower === 'compression-controlled') {
        return '<span class="status-warn">' + escapeHtml(str) + '</span>';
    }

    // Numbers: round display
    if (typeof val === 'number') {
        if (Math.abs(val) < 0.001 && val !== 0) {
            return val.toExponential(4);
        }
        if (Number.isInteger(val)) return String(val);
        return val.toFixed(4);
    }

    return escapeHtml(str);
}

/**
 * Escape HTML entities.
 */
function escapeHtml(str) {
    var div = document.createElement('div');
    div.textContent = str;
    return div.innerHTML;
}

/**
 * Show a brief toast notification.
 */
function showToast(msg, type) {
    var bg = type === 'error' ? '#dc2626' : type === 'warn' ? '#d97706' : '#16a34a';
    var t = document.createElement('div');
    t.className = 'toast';
    t.textContent = msg;
    t.style.background = bg;
    document.body.appendChild(t);
    requestAnimationFrame(function() { t.classList.add('show'); });
    setTimeout(function() {
        t.classList.remove('show');
        setTimeout(function() { t.remove(); }, 300);
    }, 4000);
}

/**
 * Fetch a calculation-sheet report and render it with CalcSheet.
 *
 * This is glue, deliberately thin: it does the network and hands the payload
 * straight to the renderer. Every number, verdict and QAQC comparison in the
 * result was decided on the server. Nothing here recomputes anything, and
 * nothing here should start to -- that is how a codebase ends up with two
 * implementations that disagree.
 *
 * @param {string} url        report endpoint
 * @param {object} payload    request body
 * @param {string} targetId   element id to render into
 * @param {object} opts       masthead fields: project, member, date
 * @param {function} done     optional callback(result)
 */
function mountSheet(url, payload, targetId, opts, done) {
    var target = document.getElementById(targetId);
    if (!target) { return; }
    target.innerHTML = '<div class="loading">Building calculation sheet...</div>';

    fetch(url, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
    })
    .then(function(res) { return res.json(); })
    .then(function(data) {
        if (data.status === 'error') {
            target.innerHTML = '<div class="error-msg">Error: '
                             + escapeHtml(data.message) + '</div>';
            return;
        }
        var o = opts || {};
        if (!o.date) { o.date = new Date().toISOString().slice(0, 10); }
        target.innerHTML = CalcSheet.render(data.result, o);
        if (window.MathJax && MathJax.typesetPromise) {
            MathJax.typesetPromise([target]);
        } else if (window.MathJax && MathJax.typeset) {
            MathJax.typeset([target]);
        }
        if (done) { done(data.result); }
    })
    .catch(function(err) {
        target.innerHTML = '<div class="error-msg">Request failed: '
                         + escapeHtml(err.message) + '</div>';
    });
}
