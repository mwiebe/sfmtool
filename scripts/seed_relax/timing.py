# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Where the chain's stages spent their wall time, reported to whoever asked.

The package still reads no clock into its own record: a stage's census is a
count of what it decided, and nothing here touches it.  What the clock reaches
instead is a SINK the caller installs -- one function taking ``(name,
seconds)`` -- so the run that owns the product owns the breakdown, and the
package holds no totals of its own.  With no sink installed the timer measures
nothing and every stage runs exactly as it ships.

``bundle_adjust`` is the native adjustment behind that same clock.  The chain
adjusts in four places (the lens on the bearings, the relaxation's own staged
adjustment, the held adjustments of the fill-in and the hand-over, and the late
release), so which stage is expensive does not by itself say whether the
adjustment is the expense: every one of them imports the adjustment from here,
and it accumulates under one name across all of them.
"""

from __future__ import annotations

import time
from contextlib import contextmanager

#: Where measured stages are reported, or ``None`` while nobody is asking.
_SINK = None


def install(sink):
    """Report every measured stage to ``sink(name, seconds)``.

    ``None`` stops the reporting.  Returns the sink that was installed before,
    so a caller that wants its own back can put it there."""
    global _SINK
    prev, _SINK = _SINK, sink
    return prev


@contextmanager
def stage(name):
    """Report the wall time of what runs inside, under ``name``.

    With nobody asking, not even the clock is read."""
    if _SINK is None:
        yield
        return
    t0 = time.perf_counter()
    try:
        yield
    finally:
        _SINK(name, time.perf_counter() - t0)


def bundle_adjust(*args, **kwargs):
    """The native bundle adjustment, timed under ``relax.bundle_adjust``.

    Imported by every module in the package that adjusts, in place of the
    kernel itself, so the adjustment's own wall time is one bucket wherever in
    the chain it was paid."""
    from sfmtool._sfmtool.geometry import bundle_adjust as _bundle_adjust

    with stage("relax.bundle_adjust"):
        return _bundle_adjust(*args, **kwargs)
