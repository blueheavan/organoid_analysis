"""Launch a Streamlit app on a free port, then shut it down once the browser tab closes.

Usage:
    python scripts/streamlit_run.py src/organoid_analysis/web_interface/organoid_workspace.py [--start-port 8501] [--idle-timeout 10]
"""
from __future__ import annotations

import argparse
import socket
import subprocess
import sys
import time
import webbrowser


def find_free_port(start: int, host: str = "127.0.0.1", span: int = 1000) -> int:
    for port in range(start, start + span):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
            probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                probe.bind((host, port))
                return port
            except OSError:
                continue
    raise RuntimeError(f"No free port found in [{start}, {start + span})")


def has_browser_connection(port: int) -> bool:
    """Whether any established TCP connection is currently attached to ``port``.

    Shells out to ``lsof`` scoped to one port rather than using psutil's
    system-wide connection scan, which raises ``AccessDenied`` on macOS when
    it reaches a process owned by another user.
    """
    result = subprocess.run(
        ["lsof", "-nP", f"-iTCP:{port}", "-sTCP:ESTABLISHED"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    return bool(result.stdout.strip())


def wait_until_serving(port: int, timeout: float) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=0.5):
                return
        except OSError:
            time.sleep(0.3)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", help="Streamlit entry-point script, e.g. src/organoid_analysis/web_interface/organoid_workspace.py")
    parser.add_argument("--start-port", type=int, default=8501)
    parser.add_argument("--startup-timeout", type=float, default=30.0,
                         help="Seconds to wait for the server to start accepting connections")
    parser.add_argument("--grace-period", type=float, default=20.0,
                         help="Seconds to allow before the first browser connection must appear")
    parser.add_argument("--idle-timeout", type=float, default=10.0,
                         help="Seconds with no browser connection before the server shuts down")
    args, extra = parser.parse_known_args()

    port = find_free_port(args.start_port)
    url = f"http://localhost:{port}"
    process = subprocess.Popen([
        sys.executable, "-m", "streamlit", "run", args.target,
        "--server.port", str(port),
        "--browser.gatherUsageStats", "false",
        "--server.headless", "true",
        *extra,
    ])
    print(f"Streamlit starting on {url} (pid {process.pid})")

    wait_until_serving(port, args.startup_timeout)
    webbrowser.open(url)

    last_connected = time.monotonic() + args.grace_period
    try:
        while process.poll() is None:
            time.sleep(2)
            if has_browser_connection(port):
                last_connected = time.monotonic()
            elif time.monotonic() - last_connected > args.idle_timeout:
                print("No browser tab attached; shutting the server down.")
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                break
    except KeyboardInterrupt:
        process.terminate()
    sys.exit(process.returncode or 0)


if __name__ == "__main__":
    main()
