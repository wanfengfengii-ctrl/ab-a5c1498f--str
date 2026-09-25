#!/usr/bin/env python3
"""HTTP 冒烟测试：对运行中的服务执行端到端检查。

通过环境变量 BASE_URL 指定服务地址（默认 http://localhost:8000，
compose 的 verify 服务中设为 http://web:8000）。
全部检查通过时以退出码 0 结束，否则以退出码 1 报告失败。
"""

import json
import os
import sys
import time
import urllib.error
import urllib.request

BASE = os.environ.get("BASE_URL", "http://localhost:8000").rstrip("/")
TIMEOUT = 30

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))


def request(method, path, payload=None):
    """返回 (status, body)；JSON 响应解析为对象，其余返回原文。"""
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(
        BASE + path,
        data=data,
        method=method,
        headers={"Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            raw = resp.read().decode()
            return resp.status, _maybe_json(raw)
    except urllib.error.HTTPError as e:
        return e.code, _maybe_json(e.read().decode())


def _maybe_json(text):
    try:
        return json.loads(text)
    except (ValueError, TypeError):
        return text


UNIQUE_PAYLOAD = {
    "candidates": [
        {"id": "C1", "genotypes": {"L1": ["12", "13"], "L2": ["7", "9"], "L3": ["16", "17"]}},
        {"id": "C2", "genotypes": {"L1": ["11", "12"], "L2": ["8", "8"], "L3": ["15", "18"]}},
        {"id": "C3", "genotypes": {"L1": ["14", "14"], "L2": ["9", "10"], "L3": ["17", "18"]}},
        {"id": "C4", "genotypes": {"L1": ["13", "15"], "L2": ["7", "10"], "L3": ["16", "16"]}},
        {"id": "C5", "genotypes": {"L1": ["11", "15"], "L2": ["8", "11"], "L3": ["15", "15"]}},
    ],
    "loci": ["L1", "L2", "L3"],
    "observed": {
        "L1": {"12": 2, "13": 2, "14": 2},
        "L2": {"7": 2, "9": 3, "10": 1},
        "L3": {"16": 2, "17": 3, "18": 1},
    },
    "tolerance": 0,
}

AMBIGUOUS_PAYLOAD = {
    "candidates": [
        {"id": "P1", "genotypes": {"L1": ["9", "10"], "L2": ["11", "12"], "L3": ["13", "14"]}},
        {"id": "P2", "genotypes": {"L1": ["9", "10"], "L2": ["11", "12"], "L3": ["13", "14"]}},
        {"id": "P3", "genotypes": {"L1": ["8", "8"], "L2": ["8", "8"], "L3": ["8", "8"]}},
        {"id": "P4", "genotypes": {"L1": ["7", "7"], "L2": ["7", "7"], "L3": ["7", "7"]}},
    ],
    "loci": ["L1", "L2", "L3"],
    "observed": {
        "L1": {"9": 2, "10": 2},
        "L2": {"11": 2, "12": 2},
        "L3": {"13": 2, "14": 2},
    },
    "tolerance": 0,
}

NO_SOLUTION_PAYLOAD = {
    **UNIQUE_PAYLOAD,
    "observed": {**UNIQUE_PAYLOAD["observed"], "L1": {"99": 6}},
}

INVALID_PAYLOAD = {
    **UNIQUE_PAYLOAD,
    "observed": {
        **UNIQUE_PAYLOAD["observed"],
        "L2": {"7": 5, "9": 3, "10": 1},  # 总深度 9 ≠ 6
    },
}


def main():
    # 等待服务就绪（最多 60 秒）
    deadline = time.time() + 60
    ready = False
    while time.time() < deadline:
        try:
            status, body = request("GET", "/api/health")
            if status == 200 and isinstance(body, dict) and body.get("status") == "ok":
                ready = True
                break
        except urllib.error.URLError:
            pass
        time.sleep(1)
    check("GET /api/health 就绪", ready)
    if not ready:
        return report()

    status, body = request("GET", "/")
    check("GET / 返回页面", status == 200 and isinstance(body, str) and "裁决" in body)

    status, body = request("POST", "/api/mixture-adjudications", UNIQUE_PAYLOAD)
    check("唯一解场景 HTTP 200", status == 200, f"实际 {status}: {body}")
    if status == 200:
        check(
            "唯一解场景判定为 unique",
            body.get("status") == "unique",
            f"实际 {body.get('status')}",
        )
        sol = body.get("solution") or {}
        check(
            "唯一解场景供样人为 C1×4 + C3×2",
            sol.get("contributors")
            == [{"candidateId": "C1", "copies": 4}, {"candidateId": "C3", "copies": 2}],
            f"实际 {sol.get('contributors')}",
        )
        check("唯一解场景残差和为 0", body.get("objective", {}).get("residualSum") == 0)

    status, body = request("POST", "/api/mixture-adjudications", AMBIGUOUS_PAYLOAD)
    check("歧义场景 HTTP 200", status == 200, f"实际 {status}: {body}")
    if status == 200:
        check(
            "歧义场景判定为 ambiguous",
            body.get("status") == "ambiguous",
            f"实际 {body.get('status')}",
        )
        witnesses = body.get("witnesses") or []
        check(
            "歧义场景给出两份不同见证",
            len(witnesses) == 2
            and witnesses[0].get("contributors") != witnesses[1].get("contributors"),
            f"实际 {witnesses}",
        )
        check("歧义场景不伪装唯一解", body.get("solution") is None)

    status, body = request("POST", "/api/mixture-adjudications", NO_SOLUTION_PAYLOAD)
    check("无解场景 HTTP 200", status == 200, f"实际 {status}: {body}")
    if status == 200:
        check(
            "无解场景判定为 no_solution",
            body.get("status") == "no_solution",
            f"实际 {body.get('status')}",
        )

    status, body = request("POST", "/api/mixture-adjudications", INVALID_PAYLOAD)
    check("非法输入返回 422", status == 422, f"实际 {status}: {body}")

    return report()


def report():
    failed = [c for c in CHECKS if not c[1]]
    for name, ok, detail in CHECKS:
        line = f"{'PASS' if ok else 'FAIL'}  {name}"
        if not ok and detail:
            line += f"  -- {detail}"
        print(line)
    print(f"\n冒烟结果：{len(CHECKS) - len(failed)}/{len(CHECKS)} 项通过")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
