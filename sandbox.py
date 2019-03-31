from pathlib import Path

from sketchy.evaluation import Evaluator

ev = Evaluator()
ev.plot_bootstraps(
    bootstrap_data=Path('tests/run1.tab'),
    confidence=0.95,
    display=True
)

