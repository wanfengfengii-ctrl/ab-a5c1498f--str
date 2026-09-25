#!/bin/sh
# 一次性验证服务：先跑代码测试，再对运行中的 API 做 HTTP 冒烟，最后以退出码报告结果。
set -eu

cd /srv/app

echo "=== [1/2] 代码测试（pytest） ==="
python -m pytest tests -q

echo "=== [2/2] HTTP 冒烟（${APP_URL:-http://app:8000}） ==="
python verify/smoke.py

echo "=== VERIFY OK ==="
