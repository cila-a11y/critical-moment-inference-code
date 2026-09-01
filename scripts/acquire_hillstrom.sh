#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DESTINATION="${1:-$ROOT/data/hillstrom_2008_raw.csv}"
LOCAL_SOURCE="${2:-}"
EXPECTED_BYTES="3964977"
EXPECTED_SHA256="0e5893329d8b93cefecc571777672028290ab69865718020c78c7284f291aece"
HTTPS_URL="https://raw.githubusercontent.com/AdityaDabrase/ab-testing-email-marketing/5866d7e8b50f46239c80d0ffa543fda501939ecc/data/raw/hillstrom.csv"

command -v sha256sum >/dev/null 2>&1 || {
  echo "GNU sha256sum is required." >&2
  exit 69
}

verify_file() {
  local path="$1"
  test "$(wc -c < "$path" | tr -d '[:space:]')" = "$EXPECTED_BYTES" &&
    test "$(sha256sum "$path" | awk '{print $1}')" = "$EXPECTED_SHA256"
}

if test -f "$DESTINATION"; then
  verify_file "$DESTINATION" || {
    echo "The existing Hillstrom file fails the locked integrity check." >&2
    exit 65
  }
  echo "HILLSTROM DATA ACQUISITION: PASSED"
  exit 0
fi

mkdir -p "$(dirname "$DESTINATION")"
TEMPORARY="$(mktemp "$(dirname "$DESTINATION")/.hillstrom.XXXXXX")"
trap 'rm -f "$TEMPORARY"' EXIT

if test -n "$LOCAL_SOURCE"; then
  test -f "$LOCAL_SOURCE" || {
    echo "Local source not found: $LOCAL_SOURCE" >&2
    exit 66
  }
  cp "$LOCAL_SOURCE" "$TEMPORARY"
else
  command -v curl >/dev/null 2>&1 || {
    echo "curl is required when no local source is supplied." >&2
    exit 69
  }
  curl --fail --location --silent --show-error \
    --proto '=https' --proto-redir '=https' --tlsv1.2 \
    --connect-timeout 20 --max-time 180 --max-filesize 3964977 \
    "$HTTPS_URL" --output "$TEMPORARY"
fi

verify_file "$TEMPORARY" || {
  echo "The downloaded Hillstrom file failed the locked integrity check." >&2
  exit 65
}
chmod 0644 "$TEMPORARY"
mv "$TEMPORARY" "$DESTINATION"
trap - EXIT
echo "HILLSTROM DATA ACQUISITION: PASSED"
