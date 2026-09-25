"""HTTP 冒烟测试：对运行中的服务执行健康检查与端到端裁决验证。

通过环境变量 APP_URL 指定服务地址（默认 http://app:8000，对应 compose 网络）。
全部检查通过以退出码 0 结束，否则退出码 1。
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.error
import urllib.request

BASE = os.environ.get("APP_URL", "http://app:8000").rstrip("/")

FAILURES: list[str] = []


def request(method: str, path: str, payload: dict | None = None):
    data = None
    headers = {}
    if payload is not None:
        data = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(BASE + path, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            raw = resp.read().decode()
            return resp.status, _maybe_json(raw)
    except urllib.error.HTTPError as exc:
        raw = exc.read().decode()
        return exc.code, _maybe_json(raw)


def _maybe_json(raw: str):
    try:
        return json.loads(raw)
    except ValueError:
        return raw


def check(name: str, condition: bool, detail: str = "") -> None:
    mark = "PASS" if condition else "FAIL"
    print(f"[{mark}] {name}" + (f" — {detail}" if detail and not condition else ""))
    if not condition:
        FAILURES.append(name)


def wait_ready(timeout: int = 60) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            status, _ = request("GET", "/api/health")
            if status == 200:
                return True
        except OSError:
            pass
        time.sleep(1)
    return False


def main() -> int:
    print(f"smoke target: {BASE}")
    check("服务就绪（GET /api/health）", wait_ready(), "服务在 60 秒内未就绪")

    status, body = request("GET", "/api/health")
    check("健康检查返回 ok", status == 200 and body == {"status": "ok"}, f"{status} {body}")

    status, body = request("GET", "/")
    check(
        "前端页面可访问",
        status == 200 and isinstance(body, str) and "裁决" in body,
        f"{status}",
    )

    unique = {
        "tolerance": 0,
        "candidates": [
            {"id": "C1", "genotypes": {"L1": ["a", "b"], "L2": ["x", "y"], "L3": ["p", "q"]}},
            {"id": "C2", "genotypes": {"L1": ["a", "c"], "L2": ["x", "x"], "L3": ["p", "p"]}},
            {"id": "C3", "genotypes": {"L1": ["b", "c"], "L2": ["y", "z"], "L3": ["q", "r"]}},
            {"id": "C4", "genotypes": {"L1": ["c", "c"], "L2": ["z", "z"], "L3": ["r", "r"]}},
        ],
        "loci": [
            {"name": "L1", "observed": {"a": 1, "b": 2, "c": 1}},
            {"name": "L2", "observed": {"x": 1, "y": 2, "z": 1}},
            {"name": "L3", "observed": {"p": 1, "q": 2, "r": 1}},
        ],
    }
    status, body = request("POST", "/api/mixture-adjudications", unique)
    ok = (
        status == 200
        and isinstance(body, dict)
        and body.get("status") == "unique"
        and body.get("objective", {}).get("contributors") == 2
        and body.get("solution", {}).get("contributors")
        == [{"id": "C1", "copies": 2}, {"id": "C3", "copies": 2}]
    )
    check("唯一解裁决（POST /api/mixture-adjudications）", ok, f"{status} {body}")

    ambiguous = {
        "tolerance": 0,
        "candidates": [
            {"id": "C1", "genotypes": {"L1": ["a", "b"], "L2": ["x", "y"], "L3": ["p", "q"]}},
            {"id": "C2", "genotypes": {"L1": ["a", "b"], "L2": ["x", "y"], "L3": ["p", "q"]}},
            {"id": "C3", "genotypes": {"L1": ["a", "c"], "L2": ["x", "z"], "L3": ["p", "r"]}},
            {"id": "C4", "genotypes": {"L1": ["c", "c"], "L2": ["z", "z"], "L3": ["r", "r"]}},
        ],
        "loci": [
            {"name": "L1", "observed": {"a": 2, "b": 2}},
            {"name": "L2", "observed": {"x": 2, "y": 2}},
            {"name": "L3", "observed": {"p": 2, "q": 2}},
        ],
    }
    status, body = request("POST", "/api/mixture-adjudications", ambiguous)
    witnesses = isinstance(body, dict) and body.get("witnesses")
    ok = (
        status == 200
        and body.get("status") == "ambiguous"
        and isinstance(witnesses, list)
        and len(witnesses) >= 2
        and witnesses[0]["contributors"] != witnesses[1]["contributors"]
    )
    check("歧义裁决出示两份不同见证", ok, f"{status} {body}")

    status, body = request(
        "POST", "/api/mixture-adjudications", {"tolerance": 0, "candidates": [], "loci": []}
    )
    check("非法输入返回 422", status == 422, f"{status} {body}")

    if FAILURES:
        print(f"\nSMOKE FAILED: {len(FAILURES)} 项未通过")
        return 1
    print("\nSMOKE OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
