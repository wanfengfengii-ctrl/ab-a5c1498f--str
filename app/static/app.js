'use strict';

const $ = (id) => document.getElementById(id);
let locusNames = [];

/* 内置示例：真实解为 C1×4 + C3×2（总深度 6，误差上限 0），有唯一解。 */
const EXAMPLE = {
  numCandidates: 5,
  numLoci: 3,
  tolerance: 0,
  candidates: [
    { id: 'C1', g: { L1: ['12', '13'], L2: ['7', '9'], L3: ['16', '17'] } },
    { id: 'C2', g: { L1: ['11', '12'], L2: ['8', '8'], L3: ['15', '18'] } },
    { id: 'C3', g: { L1: ['14', '14'], L2: ['9', '10'], L3: ['17', '18'] } },
    { id: 'C4', g: { L1: ['13', '15'], L2: ['7', '10'], L3: ['16', '16'] } },
    { id: 'C5', g: { L1: ['11', '15'], L2: ['8', '11'], L3: ['15', '15'] } },
  ],
  observed: {
    L1: { '12': 2, '13': 2, '14': 2 },
    L2: { '7': 2, '9': 3, '10': 1 },
    L3: { '16': 2, '17': 3, '18': 1 },
  },
};

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
  }[c]));
}

function fmt(x) {
  return Number.isInteger(x) ? String(x) : Number(x).toFixed(1);
}

/* ---------------- 表单构建 ---------------- */

function buildForm() {
  const nc = Number($('num-candidates').value);
  const nl = Number($('num-loci').value);
  locusNames = Array.from({ length: nl }, (_, i) => 'L' + (i + 1));

  let html = '<table class="grid"><thead><tr><th rowspan="2">编号</th>';
  for (const l of locusNames) html += `<th colspan="2">${escapeHtml(l)}</th>`;
  html += '</tr><tr>';
  for (const l of locusNames) html += '<th>等位基因 1</th><th>等位基因 2</th>';
  html += '</tr></thead><tbody id="candidate-tbody">';
  for (let i = 0; i < nc; i++) {
    html += `<tr><td><input class="cid" value="C${i + 1}" maxlength="32"></td>`;
    for (const l of locusNames) {
      html += `<td><input class="allele" data-locus="${l}" maxlength="24"></td>`;
      html += `<td><input class="allele" data-locus="${l}" maxlength="24"></td>`;
    }
    html += '</tr>';
  }
  html += '</tbody></table>';
  $('candidate-table').innerHTML = html;

  const grid = $('observed-grid');
  grid.innerHTML = '';
  for (const l of locusNames) grid.appendChild(buildLocusCard(l, {}));
  updateDepthStatus();
  $('result-panel').hidden = true;
}

function buildLocusCard(locus, table) {
  const card = document.createElement('div');
  card.className = 'locus-card';
  card.id = `obs-${locus}`;
  card.dataset.locus = locus;
  let rows = '';
  const entries = Object.entries(table);
  const list = entries.length ? entries : [['', ''], ['', ''], ['', '']];
  for (const [allele, count] of list) rows += observedRow(allele, count);
  card.innerHTML = `
    <h3>${escapeHtml(locus)} <span class="depth">总深度 0</span></h3>
    <table>
      <thead><tr><th>等位基因</th><th>观测拷贝数</th><th></th></tr></thead>
      <tbody>${rows}</tbody>
    </table>
    <button type="button" class="add-allele">+ 添加等位基因</button>`;
  return card;
}

function observedRow(allele, count) {
  return `<tr>
    <td><input class="obs-allele" value="${escapeHtml(allele)}" maxlength="24"></td>
    <td><input class="obs-count" type="number" min="0" step="1" value="${count === '' ? '' : Number(count)}"></td>
    <td><button type="button" class="del" title="删除">×</button></td>
  </tr>`;
}

function onObservedClick(ev) {
  if (ev.target.classList.contains('add-allele')) {
    const tbody = ev.target.closest('.locus-card').querySelector('tbody');
    tbody.insertAdjacentHTML('beforeend', observedRow('', ''));
  } else if (ev.target.classList.contains('del')) {
    ev.target.closest('tr').remove();
    updateDepthStatus();
  }
}

function locusTotal(locus) {
  let total = 0;
  document.querySelectorAll(`#obs-${locus} .obs-count`).forEach((inp) => {
    const n = Number(inp.value);
    if (inp.value !== '' && Number.isFinite(n)) total += n;
  });
  return total;
}

function updateDepthStatus() {
  const totals = locusNames.map((l) => {
    const t = locusTotal(l);
    const span = document.querySelector(`#obs-${l} .depth`);
    if (span) span.textContent = `总深度 ${t}`;
    return t;
  });
  const el = $('depth-status');
  const same = totals.every((t) => t === totals[0]);
  el.textContent = `各位点观测总深度：${locusNames.map((l, i) => `${l}=${totals[i]}`).join('，')} ` +
    (same ? '✓ 一致' : '✗ 不一致（裁决要求各位点总深度相同）');
  el.className = 'depth-status ' + (same ? 'ok' : 'bad');
}

/* ---------------- 示例 ---------------- */

function loadExample() {
  $('num-candidates').value = String(EXAMPLE.numCandidates);
  $('num-loci').value = String(EXAMPLE.numLoci);
  $('tolerance').value = String(EXAMPLE.tolerance);
  buildForm();
  const rows = document.querySelectorAll('#candidate-tbody tr');
  EXAMPLE.candidates.forEach((c, i) => {
    const tr = rows[i];
    tr.querySelector('.cid').value = c.id;
    locusNames.forEach((l) => {
      const inputs = tr.querySelectorAll(`input[data-locus="${l}"]`);
      inputs[0].value = c.g[l][0];
      inputs[1].value = c.g[l][1];
    });
  });
  const grid = $('observed-grid');
  grid.innerHTML = '';
  for (const l of locusNames) grid.appendChild(buildLocusCard(l, EXAMPLE.observed[l] || {}));
  updateDepthStatus();
}

/* ---------------- 提交与渲染 ---------------- */

function collectPayload() {
  const candidates = [];
  document.querySelectorAll('#candidate-tbody tr').forEach((tr) => {
    const genotypes = {};
    locusNames.forEach((l) => {
      const inputs = tr.querySelectorAll(`input[data-locus="${l}"]`);
      genotypes[l] = [inputs[0].value.trim(), inputs[1].value.trim()];
    });
    candidates.push({ id: tr.querySelector('.cid').value.trim(), genotypes });
  });
  const observed = {};
  locusNames.forEach((l) => {
    const table = {};
    document.querySelectorAll(`#obs-${l} tbody tr`).forEach((tr) => {
      const allele = tr.querySelector('.obs-allele').value.trim();
      const raw = tr.querySelector('.obs-count').value;
      if (allele && raw !== '') table[allele] = Number(raw);
    });
    observed[l] = table;
  });
  return {
    candidates,
    loci: [...locusNames],
    observed,
    tolerance: Number($('tolerance').value),
  };
}

async function adjudicate() {
  const payload = collectPayload();
  $('result-panel').hidden = false;
  $('result').innerHTML = '<p class="muted">裁决计算中…</p>';
  let resp, data;
  try {
    resp = await fetch('/api/mixture-adjudications', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    data = await resp.json();
  } catch (e) {
    $('result').innerHTML = `<div class="banner error">请求失败：${escapeHtml(e.message)}</div>`;
    return;
  }
  if (!resp.ok) renderError(resp.status, data);
  else renderResult(data);
  $('result-panel').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
}

function renderError(status, data) {
  let items = [];
  if (data && Array.isArray(data.detail)) {
    items = data.detail.map((d) => (typeof d === 'string' ? d : d.msg || JSON.stringify(d)));
  } else if (data && data.detail) {
    items = [String(data.detail)];
  } else {
    items = [JSON.stringify(data)];
  }
  $('result').innerHTML =
    `<div class="banner error">输入非法（HTTP ${status}），请修正后重新裁决：</div>` +
    `<ul class="error-list">${items.map((m) => `<li>${escapeHtml(m)}</li>`).join('')}</ul>`;
}

function renderResult(data) {
  const obj = data.objective;
  const meta =
    `<p class="objective">观测总深度 <strong>${data.totalDepth}</strong>　` +
    `误差上限 <strong>${data.tolerance}</strong>` +
    (obj
      ? `　最优供样人数 <strong>${obj.contributors}</strong>　绝对残差和 <strong>${fmt(obj.residualSum)}</strong>` +
        `　并列最优份数向量 <strong>${data.optimalVectorCount}</strong> 个`
      : '') +
    '</p>';

  if (data.status === 'unique') {
    $('result').innerHTML =
      '<div class="banner ok">唯一解：存在且仅存在一种最优供样人组合与份数向量。</div>' +
      meta + renderSolution(data.solution, '最优解');
  } else if (data.status === 'ambiguous') {
    $('result').innerHTML =
      `<div class="banner warn">歧义：同一最优目标值存在 ${data.optimalVectorCount} 种不同份数向量，` +
      `以下列出两份不同见证，不得按唯一解采信。</div>` +
      meta +
      '<div class="witness-row">' +
      renderSolution(data.witnesses[0], '见证 A') +
      renderSolution(data.witnesses[1], '见证 B') +
      '</div>';
  } else {
    $('result').innerHTML =
      '<div class="banner info">无解：在 1–3 名供样人、份数和等于观测总深度且全部残差不超过误差上限的' +
      '联合空间中，不存在任何可行组合。</div>' + meta;
  }
}

function renderSolution(sol, title) {
  let html = `<div class="solution"><h3>${escapeHtml(title)}` +
    `<span class="tag">绝对残差和 ${fmt(sol.residualSum)}</span></h3>`;
  html += '<table><thead><tr><th>候选人编号</th><th>模板份数</th></tr></thead><tbody>';
  for (const c of sol.contributors) {
    html += `<tr><td>${escapeHtml(c.candidateId)}</td><td class="num">${c.copies}</td></tr>`;
  }
  html += '</tbody></table>';
  for (const locus of sol.loci) {
    html += `<h4>位点 ${escapeHtml(locus.locus)}</h4>` +
      '<table><thead><tr><th>等位基因</th><th>观测</th><th>预测</th><th>残差</th></tr></thead><tbody>';
    for (const a of locus.alleles) {
      const cls = a.residual !== 0 ? ' class="nonzero"' : '';
      html += `<tr${cls}><td>${escapeHtml(a.allele)}</td>` +
        `<td class="num">${fmt(a.observed)}</td>` +
        `<td class="num">${fmt(a.predicted)}</td>` +
        `<td class="num">${fmt(a.residual)}</td></tr>`;
    }
    html += '</tbody></table>';
  }
  return html + '</div>';
}

/* ---------------- 初始化 ---------------- */

document.addEventListener('DOMContentLoaded', () => {
  $('btn-build').addEventListener('click', buildForm);
  $('btn-example').addEventListener('click', loadExample);
  $('btn-adjudicate').addEventListener('click', adjudicate);
  $('observed-grid').addEventListener('click', onObservedClick);
  $('observed-grid').addEventListener('input', updateDepthStatus);
  buildForm();
});
