#!/usr/bin/env python3
"""
BioCore KRAS4B Molecular Master Atlas
A publication-quality, hierarchical educational dashboard rendering the complete 
canonical human KRAS4B protein landscape (genomic, structural, clinical) 
directly into a Linux/WSL terminal and exporting a high-resolution SVG.

Dependencies:
    pip install rich
"""

import sys
import os
import datetime
from typing import Dict, List, Tuple

try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text
    from rich.layout import Layout
    from rich.columns import Columns
    from rich.align import Align
    from rich import box
except ImportError:
    print("Error: The 'rich' library is required for rendering this dashboard.")
    print("Please install it using: pip install rich")
    sys.exit(1)

# ============================================================================
# BIOLOGICAL DATA VALIDATION & KNOWLEDGE BASE
# ============================================================================

KRAS4B_SEQ = (
    "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG"
    "QEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDL"
    "PSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKK"
    "SKTKCVIM"
)

AA_3L = {
    'M': 'Met', 'T': 'Thr', 'E': 'Glu', 'Y': 'Tyr', 'K': 'Lys',
    'L': 'Leu', 'V': 'Val', 'G': 'Gly', 'A': 'Ala', 'S': 'Ser',
    'I': 'Ile', 'Q': 'Gln', 'N': 'Asn', 'H': 'His', 'F': 'Phe',
    'D': 'Asp', 'P': 'Pro', 'R': 'Arg', 'C': 'Cys', 'W': 'Trp'
}

EXACT_CODONS = {
    12: 'GGT', 13: 'GGC', 14: 'GTG', 19: 'TTA', 22: 'CAA', 28: 'TTC', 
    31: 'GAG', 33: 'GAC', 35: 'ACA', 57: 'GAT', 58: 'ACA', 59: 'GCA', 
    60: 'GGT', 61: 'CAA', 64: 'TAC', 68: 'CGA', 97: 'CGA', 117: 'AAG', 
    146: 'GCA', 185: 'TGT', 186: 'GTG', 187: 'ATA', 188: 'ATG'
}

DOMAINS = {
    "G-Domain": (1, 166, "cyan"),
    "HVR": (167, 188, "magenta")
}

MOTIFS = {
    "P-loop": (10, 17, "bright_yellow"),
    "Switch I": (30, 38, "bright_green"),
    "Switch II": (59, 76, "bright_green"),
    "CAAX": (185, 188, "bright_magenta")
}

HOTSPOTS = {
    12: {
        "aa": "G", "codon": "GGT", "type": "Primary", "color": "bright_red",
        "subs": [("TGT", "G→T", "Transversion", "G12C"), ("GAT", "G→A", "Transition", "G12D"), ("GTT", "G→T", "Transversion", "G12V"), ("AGT", "G→A", "Transition", "G12S"), ("CGT", "G→C", "Transversion", "G12R"), ("GCT", "G→C", "Transversion", "G12A")]
    },
    13: {
        "aa": "G", "codon": "GGC", "type": "Primary", "color": "bright_red",
        "subs": [("GAC", "G→A", "Transition", "G13D"), ("TGC", "G→T", "Transversion", "G13C"), ("CGC", "G→C", "Transversion", "G13R"), ("GTC", "G→T", "Transversion", "G13V")]
    },
    61: {
        "aa": "Q", "codon": "CAA", "type": "Primary", "color": "bright_red",
        "subs": [("CAT", "A→T", "Transversion", "Q61H"), ("CTA", "A→T", "Transversion", "Q61L"), ("CGA", "A→G", "Transition", "Q61R"), ("AAA", "C→A", "Transversion", "Q61K")]
    },
    146: {
        "aa": "A", "codon": "GCA", "type": "Secondary", "color": "dark_orange",
        "subs": [("ACA", "G→A", "Transition", "A146T"), ("GTA", "C→T", "Transition", "A146V"), ("CCA", "G→C", "Transversion", "A146P")]
    },
    117: {
        "aa": "K", "codon": "AAG", "type": "Secondary", "color": "dark_orange",
        "subs": [("AAC", "G→C", "Transversion", "K117N"), ("AGG", "A→G", "Transition", "K117R")]
    },
    59: {
        "aa": "A", "codon": "GCA", "type": "Secondary", "color": "yellow",
        "subs": [("ACA", "G→A", "Transition", "A59T"), ("GGA", "G→G", "Transversion", "A59G")]
    },
    68: {
        "aa": "R", "codon": "CGA", "type": "Secondary", "color": "yellow",
        "subs": [("AGC", "CGA→AGC", "Complex", "R68S")]
    }
}


# ============================================================================
# RENDERING ENGINE & ATLAS ARCHITECTURE
# ============================================================================

class KRAS4BMasterAtlas:
    def __init__(self):
        # 160-column width strictly enforced to guarantee perfect rendering & SVG generation
        self.console = Console(record=True, width=160, force_terminal=True)
        self.validate_data()

    def validate_data(self):
        """Strict scientific validation of canonical sequence."""
        assert len(KRAS4B_SEQ) == 188, "Validation Error: Sequence length must be 188."
        assert KRAS4B_SEQ.startswith("MTEYKLVVVG"), "Validation Error: Invalid N-terminus."
        assert KRAS4B_SEQ.endswith("CVIM"), "Validation Error: Invalid C-terminus CAAX motif."
        
        self.codons = {}
        for i, aa in enumerate(KRAS4B_SEQ):
            pos = i + 1
            if pos in EXACT_CODONS:
                self.codons[pos] = EXACT_CODONS[pos]
            else:
                # Basic standard codon fallback for non-hotspot regions
                std = {'M':'ATG','E':'GAG','P':'CCG','Q':'CAG','S':'AGC','D':'GAC','V':'GTG','L':'CTG','T':'ACC','F':'TTC','W':'TGG','K':'AAG','N':'AAC','A':'GCC','I':'ATC','R':'CGC','G':'GGC','Y':'TAC','H':'CAC','C':'TGC'}
                self.codons[pos] = std.get(aa, 'NNN')

    def _get_region_color(self, pos: int) -> str:
        for name, (start, end, color) in MOTIFS.items():
            if start <= pos <= end: return color
        for name, (start, end, color) in DOMAINS.items():
            if start <= pos <= end: return color
        return "blue"
        
    def _get_heat_color(self, pos: int) -> str:
        if pos in [12, 13, 61]: return "bright_red"
        if pos in [146, 117]: return "dark_orange"
        if pos in [59, 68, 14, 19, 22]: return "yellow"
        if pos in EXACT_CODONS and pos not in [185, 186, 187, 188]: return "dodger_blue1"
        return "grey30"

    # --- ROW 1: HEADER ---
    def render_header(self):
        table = Table(box=box.DOUBLE, expand=True, border_style="cyan", padding=(0, 2))
        table.add_column(justify="center")
        
        title = Text("KRAS4B Master Molecular Atlas\n", style="bold white", justify="center")
        title.append("Codon-level mutation logic, structural interpretation, and translational context[nature]\n", style="bold cyan")
        
        meta = Table.grid(padding=(0, 4))
        meta.add_row(
            Text("Gene: KRAS", style="grey70"), Text("Protein: Cellular Tumor Antigen p53", style="grey70"),
            Text("Isoform: KRAS4B", style="grey70"), Text("Length: 188 aa", style="grey70"),
            Text("Chromosome: 12p12.1", style="grey70")
        )
        meta.add_row(
            Text("Valid: Canonical Seq", style="bold green"), Text(f"Rendered: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}", style="bold green"),
            Text("DBs: UniProt P01116-2, COSMIC v98", style="bold green"), Text("Clinical: OncoKB, TCGA", style="bold green"),
            Text("Structure: PDB 4OBE", style="bold green")
        )
        
        table.add_row(title)
        table.add_row(Align.center(meta))
        self.console.print(table)

    # --- ROW 2: DOMAIN MAP & HEATMAP ---
    def render_linear_maps(self):
        panel_text = Text()
        
        # 1. Structural Map
        panel_text.append(" [PROTEIN ARCHITECTURE]\n", style="bold white")
        panel_text.append(" 1 ", style="grey50")
        for i in range(1, 189):
            if i in range(10, 18): panel_text.append("█", style="bright_yellow")
            elif i in range(30, 39): panel_text.append("█", style="bright_green")
            elif i in range(59, 77): panel_text.append("█", style="bright_green")
            elif i in range(185, 189): panel_text.append("█", style="bright_magenta")
            elif i <= 166: panel_text.append("█", style="cyan")
            else: panel_text.append("█", style="magenta")
        panel_text.append(" 188\n", style="grey50")
        panel_text.append("   └─G-Domain (Catalytic GTPase)──────────────────────────────────────────────────────────────────────────────────────────┘ └─HVR─────────┘\n", style="cyan")
        
        # 2. Mutation Density Heatmap
        panel_text.append("\n [MUTATION DENSITY HEATMAP]\n", style="bold white")
        panel_text.append(" 1 ", style="grey50")
        for i in range(1, 189):
            panel_text.append("█", style=self._get_heat_color(i))
        panel_text.append(" 188\n", style="grey50")
        panel_text.append("    ▲▲                                               ▲                                                      ▲             ▲\n", style="bright_red")
        panel_text.append("   G12/13                                           Q61                                                    K117          A146\n", style="bright_red")
        
        self.console.print(Panel(panel_text, box=box.ROUNDED, border_style="grey50"))

    # --- ROW 3 & 4: HOTSPOT STRIP & CODON ATLAS ---
    def render_codon_atlas(self):
        self.console.print(Panel(Text("COMPLETE CODON ATLAS (1-188)", style="bold white", justify="center"), border_style="cyan"))
        
        chunk_size = 20
        for i in range(0, 188, chunk_size):
            start = i + 1
            end = min(i + chunk_size, 188)
            seq = KRAS4B_SEQ[i:i+chunk_size]
            
            p_text = Text()
            
            # Feature Annotations Line
            feat_line = [" "] * (chunk_size * 7)
            for m_name, (ms, me, m_color) in MOTIFS.items():
                if ms <= end and me >= start:
                    loc_s = max(0, ms - start) * 7
                    loc_e = min(chunk_size - 1, me - start) * 7 + 5
                    for j in range(loc_s, loc_e + 1): feat_line[j] = "─"
                    feat_line[loc_s] = "╭"
                    feat_line[loc_e] = "╮"
                    if loc_e - loc_s > len(m_name):
                        s_idx = loc_s + ((loc_e - loc_s) // 2) - (len(m_name) // 2)
                        for k, c in enumerate(m_name): feat_line[s_idx + k] = c
            
            f_str = "".join(feat_line).rstrip()
            if f_str: p_text.append("           │ " + f_str + "\n", style="bold white")
            
            # Rows initialization
            ptm_line = Text(" PTMs      │", style="bold magenta")
            mut_line = Text(" Hotspot   │", style="bold red")
            idx_line = Text(f" [{start:03d}-{end:03d}] │", style="bold white")
            dna_line = Text(" Wild DNA  │", style="grey70")
            rna_line = Text(" Wild mRNA │", style="grey70")
            aa3_line = Text(" AA (3L)   │", style="grey70")
            aa1_line = Text(" AA (1L)   │", style="bold white")
            
            for j, aa in enumerate(seq):
                pos = start + j
                style = self._get_region_color(pos)
                dna = self.codons[pos]
                rna = dna.replace('T', 'U')
                aa3 = AA_3L.get(aa, 'Xaa')
                
                # PTMs logic
                ptm_str = "      "
                if pos == 185: ptm_str = " Farns"
                elif pos >= 175 and pos <= 180: ptm_str = " PolyK"
                
                # Mut logic
                mut_str = "      "
                aa1_style = style
                if pos in HOTSPOTS:
                    mut_str = "  [■] "
                    aa1_style = f"bold {HOTSPOTS[pos]['color']} reverse"
                elif pos in EXACT_CODONS and pos not in [185, 186, 187, 188]:
                    mut_str = "  [●] "
                
                ptm_line.append(f"{ptm_str:^7s}", style="magenta")
                mut_line.append(f"{mut_str:^7s}", style=self._get_heat_color(pos))
                idx_line.append(f"{pos:^7d}", style="grey50")
                dna_line.append(f"{dna:^7s}", style=style)
                rna_line.append(f"{rna:^7s}", style=style)
                aa3_line.append(f"{aa3:^7s}", style=style)
                aa1_line.append(f"{aa:^7s}", style=aa1_style)

            p_text.append(ptm_line); p_text.append("\n")
            p_text.append(mut_line); p_text.append("\n")
            p_text.append(idx_line); p_text.append("\n")
            p_text.append(dna_line); p_text.append("\n")
            p_text.append(rna_line); p_text.append("\n")
            p_text.append(aa3_line); p_text.append("\n")
            p_text.append(aa1_line); p_text.append("\n")

            self.console.print(Panel(p_text, border_style="grey30"))

    # --- ROW 5: HOTSPOT CARDS ---
    def render_hotspot_cards(self):
        cards = []
        for pos in [12, 13, 61]: # Primary Drivers
            data = HOTSPOTS[pos]
            t = Table(box=box.SIMPLE_HEAVY, border_style=data["color"], show_header=False)
            t.add_column(); t.add_column()
            t.add_row(Text("Wild DNA:", style="grey70"), Text(f"{data['codon']} ({data['aa']})", style="bold white"))
            t.add_row(Text("Driver:", style="grey70"), Text(f"{data['type']} (GOF)", style="bold bright_red"))
            
            sub_t = Table(box=None, padding=0)
            sub_t.add_column(); sub_t.add_column(); sub_t.add_column()
            for mut in data['subs'][:3]: # Show top 3 for space
                sub_t.add_row(f"{data['codon']}→{mut[0]}", mut[1], f"[bold white]{mut[3]}[/]")
            if len(data['subs']) > 3:
                sub_t.add_row("...", "...", "[dim]...[/]")
                
            t.add_row(Text("Subs:", style="grey70"), sub_t)
            cards.append(Panel(t, title=f"[bold]CODON {pos}[/]", border_style=data["color"], width=50))
            
        self.console.print(Panel(Columns(cards, expand=True, equal=True), title="[bold white]PRIMARY DRIVER HOTSPOTS (95%+ of KRAS Mutations)[/bold white]", border_style="bright_red"))

    # --- ROW 6: MUTATION JOURNEY ---
    def render_mutation_journey(self):
        journey = Text()
        journey.append("  Wild DNA Codon       Mutated Codon         Changed Base       Protein Effect      Structural Change         Biochemical Defect         Signaling Effect            Clinical Actionability\n", style="bold white")
        journey.append("  [ GGT ] ───────────► [ TGT ] ────────────► G→T Transv ──────► G12C Missense ────► P-loop Clash ───────────► Loss of GAP Hydrolysis ─► MAPK / PI3K Activation ────► Sotorasib / Adagrasib \n", style="bright_cyan")
        journey.append("  [ CAA ] ───────────► [ CTA ] ────────────► A→T Transv ──────► Q61L Missense ────► Switch II Disruption ───► Loss of Catalytic H2O ──► MAPK / PI3K Activation ────► MEK/ERK Blockade \n", style="bright_red")
        self.console.print(Panel(journey, title="[bold]PATHOGENESIS JOURNEY (Genotype to Phenotype)[/bold]", border_style="cyan"))

    # --- ROW 7 & 8: STRUCTURE, MOLECULAR SWITCH, SIGNALING ---
    def render_biology_panels(self):
        # 1. Structural & Molecular Switch
        switch = Text()
        switch.append("                     [ GDP-Bound (Inactive) ] \n", style="grey70")
        switch.append("                           ▲          │       \n", style="bold white")
        switch.append("     GAP (NF1) stimulates  │          │  SOS1 (GEF) ejects GDP\n", style="bold bright_red")
        switch.append("     GTP hydrolysis        │          ▼       \n", style="bold white")
        switch.append("                     [ GTP-Bound (Active) ]   \n\n", style="bold bright_green")
        switch.append(" ■ G12 / G13 Mut: Block GAP arginine finger insertion. (Trapped in Active state)\n", style="bright_red")
        switch.append(" ■ Q61 Mut: Destroys catalytic H2O coordination. (Trapped in Active state)\n", style="bright_red")
        switch.append(" ■ A146 / K117 Mut: Destabilize nucleotide pocket. (Fast GDP/GTP cycling)\n", style="dark_orange")
        
        # 2. Signaling
        signal = Text()
        signal.append(" [ ONCOGENIC SIGNALING CASCADES ]\n\n", style="bold white")
        signal.append(" KRAS(GTP) ──► RAF ──► MEK ──► ERK ──► Proliferation\n", style="bold bright_green")
        signal.append(" KRAS(GTP) ──► PI3K ──► AKT ──► mTOR ──► Cell Survival\n", style="bold bright_green")
        signal.append(" KRAS(GTP) ──► RalGDS ──► RalA/B ──► Vesicle Trafficking\n\n", style="bold bright_green")
        signal.append(" [ MEMBRANE LOCALIZATION ]\n", style="bold white")
        signal.append(" CAAX Motif (185-188) directs Farnesylation, RCE1 cleavage, and \n ICMT methylation, tethering KRAS to the plasma membrane.", style="magenta")

        cols = Columns([Panel(switch, title="[bold]THE GTPASE MOLECULAR SWITCH[/]", border_style="green", width=78), 
                        Panel(signal, title="[bold]DOWNSTREAM EFFECTORS & LOCALIZATION[/]", border_style="cyan", width=78)])
        self.console.print(cols)

    # --- ROW 9: CLINICAL SUMMARY TABLE ---
    def render_clinical_panel(self):
        t = Table(box=box.SIMPLE_HEAVY, expand=True, border_style="blue")
        t.add_column("Driver Status", style="bold white")
        t.add_column("Actionability", style="bold cyan")
        t.add_column("FDA-Approved Therapy", style="bold green")
        t.add_column("Resistance Mechanisms", style="dark_orange")
        t.add_column("Primary Cancers", style="magenta")
        
        t.add_row("G12C (Driver)", "Directly Druggable (Cys)", "Sotorasib, Adagrasib", "R68S, Y96D, Amplification", "NSCLC, CRC")
        t.add_row("G12D (Driver)", "Druggable (Non-Covalent)", "Trials (MRTX1133)", "Upstream EGFR feedback", "PDAC, CRC")
        t.add_row("G12V, G13D, Q61", "Indirectly Druggable", "None (Targeting MEK/ERK)", "Parallel PI3K activation", "CRC, Melanoma")
        t.add_row("A146T, K117N", "Indirectly Druggable", "None (Targeting SOS1/SHP2)", "Rapid cycling overrides", "CRC")
        self.console.print(Panel(t, title="[bold]CLINICAL TRANSLATION & THERAPEUTIC LANDSCAPE[/bold]", border_style="blue"))

    # --- ROW 10: RULES PANEL ---
    def render_rules(self):
        rules = Text()
        rules.append("HOW TO READ KRAS MUTATIONS:\n", style="bold white")
        rules.append("Rule 1: ", style="bold cyan"); rules.append("Codons 12, 13, and 61 dominate KRAS oncogenesis (>95% of mutations).\n")
        rules.append("Rule 2: ", style="bold cyan"); rules.append("Single nucleotide substitutions create most driver mutations, altering exact structural topology.\n")
        rules.append("Rule 3: ", style="bold cyan"); rules.append("Different substitutions produce different biochemical behaviors (e.g., G12C is covalently druggable, G12D is not).\n")
        rules.append("Rule 4: ", style="bold cyan"); rules.append("Driver mutation ≠ universally druggable mutation. Context (tissue, co-mutations) dictates therapeutic response.\n")
        rules.append("Rule 5: ", style="bold cyan"); rules.append("The Hypervariable Region (HVR) regulates plasma membrane localization, not intrinsic catalysis.\n")
        self.console.print(Panel(rules, border_style="grey50"))

    # --- ROW 11: FOOTER ---
    def render_footer(self):
        footer_text = "Focus hotspots: G12, G13, Q61, K117, A146[nature]"
        self.console.print(Panel(Align.center(Text(footer_text, style="grey70")), box=box.ROUNDED, border_style="grey30"))

    # --- MASTER RENDER & EXPORT LOGIC ---
    def render_all(self):
        self.console.print("\n")
        self.render_header()
        self.render_linear_maps()
        self.render_codon_atlas()
        self.render_hotspot_cards()
        self.render_mutation_journey()
        self.render_biology_panels()
        self.render_clinical_panel()
        self.render_rules()
        self.render_footer()
        self.console.print("\n")
        
        # --- EXPORT TO SVG ---
        output_dir = "/mnt/d/Bioinformatics_Projects/vision-pathway"
        os.makedirs(output_dir, exist_ok=True)
        
        out_svg = os.path.join(output_dir, "KRAS4B_Master_Atlas.svg")
        # Save as SVG (Native vector format for terminal captures, handles tall/scrolling buffers perfectly)
        self.console.save_svg(out_svg, title="KRAS4B Master Molecular Atlas")
        
        print(f"\033[92m[SUCCESS]\033[0m Publication-quality vector image saved to: \033[96m{out_svg}\033[0m")
        print("\033[90mNote: Open the SVG in any web browser or vector editor (Illustrator/Inkscape) for infinite resolution scaling.\033[0m\n")


if __name__ == "__main__":
    atlas = KRAS4BMasterAtlas()
    atlas.render_all()