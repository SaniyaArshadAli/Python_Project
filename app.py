import streamlit as st
import pandas as pd


# ---------------------------------------------------------
# Page configuration
# ---------------------------------------------------------

st.set_page_config(
    page_title="PWM DNA Motif Scanner",
    page_icon="🧬",
    layout="wide",
)


# ---------------------------------------------------------
# Default values
# ---------------------------------------------------------

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


# ---------------------------------------------------------
# PWM scanning logic
# ---------------------------------------------------------

def scan_pwm(dna: str, pwm: dict, threshold: float):
    """
    Scan a DNA sequence using a position-weight matrix.

    Each window has MOTIF_LEN bases. The score is calculated by
    adding the PWM value for each base at each position.
    """

    matches = []
    all_windows = []

    for i in range(len(dna) - MOTIF_LEN + 1):
        window = dna[i:i + MOTIF_LEN]
        score = 0.0

        for position, base in enumerate(window):
            score += pwm[base][position]

        rounded_score = round(score, 2)
        is_match = score >= threshold

        entry = {
            "Start": i,
            "End": i + MOTIF_LEN - 1,
            "Window": window,
            "Score": rounded_score,
            "Match": is_match,
        }

        all_windows.append(entry)

        if is_match:
            matches.append(entry)

    return matches, all_windows


# ---------------------------------------------------------
# Application title
# ---------------------------------------------------------

st.title("🧬 PWM DNA Motif Scanner")

st.write(
    "Scan a DNA sequence using an editable position-weight matrix."
)

st.divider()


# ---------------------------------------------------------
# Input form
# ---------------------------------------------------------

with st.form("pwm_scanner_form"):

    st.subheader("Scan settings")

    dna_input = st.text_input(
        "DNA sequence",
        value=DEFAULT_DNA,
        help="Enter a sequence containing only A, C, G, and T.",
    )

    threshold_input = st.number_input(
        "Score threshold",
        value=DEFAULT_THRESHOLD,
        step=0.1,
        format="%.2f",
    )

    st.subheader("Position-weight matrix")

    st.write(
        "Edit the PWM values below. Each column represents one motif position."
    )

    pwm_values = {}

    columns = st.columns([1, 1, 1, 1, 1])

    columns[0].markdown("**Base**")
    columns[1].markdown("**Position 0**")
    columns[2].markdown("**Position 1**")
    columns[3].markdown("**Position 2**")
    columns[4].markdown("**Description**")

    for base in ["A", "C", "G", "T"]:
        row_columns = st.columns([1, 1, 1, 1, 1])

        row_columns[0].markdown(f"### {base}")

        pwm_values[base] = []

        for position in range(MOTIF_LEN):
            value = row_columns[position + 1].number_input(
                label=f"{base}_{position}",
                value=float(DEFAULT_PWM[base][position]),
                step=0.1,
                format="%.2f",
                key=f"pwm_{base}_{position}",
                label_visibility="collapsed",
            )

            pwm_values[base].append(value)

        descriptions = {
            "A": "Adenine",
            "C": "Cytosine",
            "G": "Guanine",
            "T": "Thymine",
        }

        row_columns[4].write(descriptions[base])

    st.divider()

    scan_button = st.form_submit_button(
        "Scan DNA",
        type="primary",
        use_container_width=True,
    )


# ---------------------------------------------------------
# Scan after form submission
# ---------------------------------------------------------

if scan_button:

    dna = dna_input.upper().strip()
    threshold = float(threshold_input)

    if not dna:
        st.error("Please enter a DNA sequence.")

    elif not set(dna).issubset(VALID_BASES):
        invalid_bases = sorted(set(dna) - VALID_BASES)

        st.error(
            "Invalid DNA sequence. Only A, C, G, and T are allowed."
        )

        st.write(
            "Invalid character(s):",
            ", ".join(invalid_bases),
        )

    elif len(dna) < MOTIF_LEN:
        st.error(
            f"DNA sequence must be at least {MOTIF_LEN} base pairs long."
        )

    else:
        matches, all_windows = scan_pwm(
            dna=dna,
            pwm=pwm_values,
            threshold=threshold,
        )

        all_windows_df = pd.DataFrame(all_windows)
        matches_df = pd.DataFrame(matches)

        st.success(
            f"Scanned {len(all_windows)} DNA window(s) and found "
            f"{len(matches)} matching window(s)."
        )

        st.divider()

        # -------------------------------------------------
        # Summary metrics
        # -------------------------------------------------

        st.subheader("Scan summary")

        metric_columns = st.columns(4)

        metric_columns[0].metric(
            "DNA length",
            len(dna),
        )

        metric_columns[1].metric(
            "Motif length",
            MOTIF_LEN,
        )

        metric_columns[2].metric(
            "Windows scanned",
            len(all_windows),
        )

        metric_columns[3].metric(
            "Matches",
            len(matches),
        )

        st.divider()

        # -------------------------------------------------
        # Matching windows
        # -------------------------------------------------

        st.subheader("Matching windows")

        if matches:
            st.dataframe(
                matches_df,
                use_container_width=True,
                hide_index=True,
            )
        else:
            st.info(
                "No DNA windows reached the selected threshold."
            )

        st.divider()

        # -------------------------------------------------
        # All windows
        # -------------------------------------------------

        st.subheader("All scanned windows")

        st.dataframe(
            all_windows_df,
            use_container_width=True,
            hide_index=True,
        )

        st.caption(
            "Start and end positions use zero-based indexing."
        )

        st.divider()

        # -------------------------------------------------
        # PWM display
        # -------------------------------------------------

        st.subheader("PWM used for this scan")

        pwm_dataframe = pd.DataFrame(
            pwm_values,
            index=[
                f"Position {position}"
                for position in range(MOTIF_LEN)
            ],
        )

        st.dataframe(
            pwm_dataframe,
            use_container_width=True,
        )
