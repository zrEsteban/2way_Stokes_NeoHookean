#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "usage: $0 MAPS_OUTPUT COMMAND [ARG ...]" >&2
    exit 2
fi

maps_output=$1
shift

if [[ -s "$maps_output" ]]; then
    exec "$@"
fi

"$@" &
child=$!

captured=no
wait_pattern=${MAPS_WAIT_PATTERN:-}
for _ in $(seq 1 500); do
    if [[ -r "/proc/$child/maps" ]] \
        && { [[ -z "$wait_pattern" ]] || grep -q "$wait_pattern" "/proc/$child/maps"; }
    then
        cp "/proc/$child/maps" "$maps_output"
        captured=yes
        break
    fi
    sleep 0.01
done

if [[ "$captured" != yes ]]; then
    echo "could not capture /proc/$child/maps" >&2
    kill "$child" 2>/dev/null || true
    wait "$child" || true
    exit 1
fi

wait "$child"
