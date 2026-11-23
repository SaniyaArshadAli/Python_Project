from flask import Flask, render_template, request

app = Flask(__name__)

DEFAULT_PWM = {
    "A": [1.5, -2.0, 0.8],
    "C": [-0.5, 1.8, -1.0],
    "G": [-0.8, -1.2, 1.9],
    "T": [-2.0, 0.5, 1.6],
}
DEFAULT_DNA = "ACGTACGATCGATCGATCG"
DEFAULT_THRESHOLD = 2.0
MOTIF_LEN = 3
VALID_BASES = {"A", "C", "G", "T"}


def scan_pwm(dna: str, pwm: dict, threshold: float):
  matches = []
  all_windows = []

  for i in range(len(dna) - MOTIF_LEN + 1):
    window = dna[i : i + MOTIF_LEN]
    score = 0.0
    for pos, base in enumerate(window):
      score += pwm[base][pos]

    rounded_score = round(score, 2)
    is_match = score >= threshold
    entry = {
        "start": i,
        "end": i + MOTIF_LEN - 1,
        "window": window,
        "score": rounded_score,
        "is_match": is_match,
    }
    all_windows.append(entry)
    if is_match:
      matches.append(entry)

  return matches, all_windows


@app.route("/", methods=["GET", "POST"])
def index():
  dna = DEFAULT_DNA
  threshold = DEFAULT_THRESHOLD
  pwm = {base: list(weights) for base, weights in DEFAULT_PWM.items()}
  error = None
  matches = []
  all_windows = []
  scanned = False

  if request.method == "POST":
    dna = request.form.get("dna", "").upper().strip()
    try:
      threshold = float(request.form.get("threshold", DEFAULT_THRESHOLD))
    except ValueError:
      threshold = DEFAULT_THRESHOLD

    # Parse editable PWM matrix from form inputs
    for base in ["A", "C", "G", "T"]:
      for pos in range(MOTIF_LEN):
        field_name = f"pwm_{base}_{pos}"
        raw_val = request.form.get(field_name)
        try:
          pwm[base][pos] = float(raw_val)
        except (TypeError, ValueError):
          pwm[base][pos] = DEFAULT_PWM[base][pos]

    if not dna or not set(dna).issubset(VALID_BASES):
      error = "Please enter a valid DNA sequence containing only A, C, G, and T."
    elif len(dna) < MOTIF_LEN:
      error = f"DNA sequence must be at least {MOTIF_LEN} base pairs long."
    else:
      matches, all_windows = scan_pwm(dna, pwm, threshold)
      scanned = True

  return render_template(
      "index.html",
      dna=dna,
      threshold=threshold,
      pwm=pwm,
      bases=["A", "C", "G", "T"],
      positions=list(range(MOTIF_LEN)),
      error=error,
      matches=matches,
      all_windows=all_windows,
      scanned=scanned,
  )


if __name__ == "__main__":
  app.run(debug=True)