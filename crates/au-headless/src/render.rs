//! Drawing what the recorder saw.
//!
//! One self-contained HTML file with inline SVG. No dependencies, no CDN, no
//! JavaScript — the simulation core has never taken a third-party crate and a
//! debugging view is a poor reason to start. It opens in a browser and it works
//! from a file:// URL on a machine with no network.
//!
//! # This is instrumentation, not rendering
//!
//! Worth being precise about, because the project's vision draws a hard line
//! here and it would be easy to blur it.
//!
//! **Rendering** — eventually — must be derived. A creature's colour will be a
//! consequence of what it is made of and how it is built; the human view will be
//! one declared mapping among possible sensors, and it will have to say so.
//! Nothing may be authored.
//!
//! **This file is not that.** These are arbitrary colours assigned to arbitrary
//! series so a person can tell two lines apart. It draws *measurements of* the
//! world, not the world. Blue here means "the second species", chosen by me,
//! meaning nothing. Confusing the two would be the first step toward authoring
//! appearance, which is the thing we have agreed not to do.
//!
//! # What it draws, and why those
//!
//! Every panel corresponds to a mistake that actually happened and went unseen:
//!
//! * **the protocell census** — a sawtooth here is a world cycling, which is what
//!   Phase 5i step 2 did for two sessions while its division counter looked
//!   healthy. Population and cumulative births/divisions/deaths together
//!   distinguish growth from churn; neither does alone.
//! * **populations, bulk against encapsulated** — since protocells, matter lives
//!   in two places, and drawing the sum would hide the only interesting thing.
//!   Log scale, because these span orders of magnitude and a linear axis shows
//!   one line and a floor.
//! * **temperature as a band** — two of the worst bugs in this project announced
//!   themselves as absurd temperatures with nobody watching. Drawn on a signed
//!   log scale so 10¹⁸ K stays on the page instead of flattening everything else
//!   into the axis.
//! * **boundary flows** — an open world's inflow and outflow, gross. Equal and
//!   rising is a steady state; equal and flat is a sealed box; diverging is
//!   accumulation.

use au_sim::observe::Sample;

const W: f64 = 900.0;
const H: f64 = 220.0;
const PAD_L: f64 = 78.0;
const PAD_R: f64 = 150.0;
const PAD_T: f64 = 26.0;
const PAD_B: f64 = 34.0;

/// Arbitrary and meaningless — see the module docs. Chosen only to be
/// distinguishable, including for the commonest colour-vision deficiencies.
const SERIES: [&str; 8] = [
    "#1f6fb4", "#d1495b", "#2a9d5c", "#c77f00", "#7052a8", "#00868b", "#a8562c", "#5a6570",
];

pub fn render(samples: &[Sample], title: &str, species: &[String]) -> String {
    let mut o = String::new();
    o.push_str(&format!(
        r#"<!doctype html><meta charset="utf-8"><title>{}</title>
<style>
 body{{background:#fbfbfa;color:#1a1a1a;font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,sans-serif;margin:0;padding:32px 40px 64px}}
 h1{{font-size:20px;font-weight:600;margin:0 0 4px}}
 .sub{{color:#6b6b6b;font-size:13px;margin-bottom:28px}}
 .panel{{background:#fff;border:1px solid #e6e4e0;border-radius:6px;margin-bottom:20px;padding:14px 16px 6px}}
 .panel h2{{font-size:13px;font-weight:600;margin:0 0 2px;letter-spacing:.01em}}
 .panel p{{color:#6b6b6b;font-size:12px;margin:0 0 6px}}
 .note{{color:#8a8a8a;font-size:12px;font-style:italic;margin-top:24px;max-width:70ch}}
 text{{font:11px -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;fill:#6b6b6b}}
 .ax{{stroke:#e6e4e0;stroke-width:1}}
</style>
<h1>{}</h1><div class="sub">"#,
        esc(title),
        esc(title)
    ));
    if let (Some(a), Some(b)) = (samples.first(), samples.last()) {
        o.push_str(&format!(
            "{} samples · tick {} → {} · {:.4} s simulated",
            samples.len(),
            a.tick,
            b.tick,
            b.seconds
        ));
    } else {
        o.push_str("no samples");
    }
    o.push_str("</div>");

    if samples.is_empty() {
        o.push_str("<p>Nothing was recorded.</p>");
        return o;
    }

    let t: Vec<f64> = samples.iter().map(|s| s.tick as f64).collect();

    // ── Protocells ──────────────────────────────────────────────────────────
    if samples.iter().any(|s| s.protocells > 0 || s.nucleated > 0) {
        let mut lines: Vec<(String, Vec<f64>)> = Vec::new();
        lines.push(("living".into(), samples.iter().map(|s| s.protocells as f64).collect()));
        lines.push(("nucleated".into(), samples.iter().map(|s| s.nucleated as f64).collect()));
        lines.push(("divided".into(), samples.iter().map(|s| s.divided as f64).collect()));
        lines.push(("lysed".into(), samples.iter().map(|s| s.lysed as f64).collect()));
        o.push_str(&panel(
            "Protocells",
            "living population, and cumulative births, divisions and deaths. \
             A sawtooth in the population with rising counters is churn, not growth.",
            &chart(&t, &lines, false),
        ));
    }

    // ── Populations ─────────────────────────────────────────────────────────
    let n_sp = samples[0].bulk.len();
    let any_enc = samples.iter().any(|s| s.encapsulated.iter().any(|&v| v != 0));
    let mut lines: Vec<(String, Vec<f64>)> = Vec::new();
    for s in 0..n_sp {
        let name = species.get(s).cloned().unwrap_or_else(|| format!("species {}", s));
        let series: Vec<f64> = samples.iter().map(|x| x.bulk[s] as f64).collect();
        if series.iter().any(|&v| v > 0.0) {
            lines.push((format!("{} · medium", name), series));
        }
        if any_enc {
            let e: Vec<f64> = samples.iter().map(|x| x.encapsulated[s] as f64).collect();
            if e.iter().any(|&v| v > 0.0) {
                lines.push((format!("{} · in bags", name), e));
            }
        }
    }
    if !lines.is_empty() {
        o.push_str(&panel(
            "Populations",
            "log scale. Medium and encapsulated are drawn apart on purpose — their \
             sum is what conservation checks, their ratio is what tells you whether \
             compartments are doing anything.",
            &chart(&t, &lines, true),
        ));
    }

    // ── Temperature ─────────────────────────────────────────────────────────
    let temps: Vec<(String, Vec<f64>)> = vec![
        ("max".into(), samples.iter().map(|s| s.temp_max).collect()),
        ("mean".into(), samples.iter().map(|s| s.temp_mean).collect()),
        ("min".into(), samples.iter().map(|s| s.temp_min).collect()),
    ];
    let hot = samples.iter().any(|s| s.temp_max > 1.0e5 || s.temp_min < 0.0);
    o.push_str(&panel(
        "Temperature (K)",
        if hot {
            "⚠ this world left any physical range. A run that does that is telling \
             you something before any test does."
        } else {
            "derived from energy and mass exactly as physics derives it, so the two \
             layers can never disagree about how hot anything is."
        },
        &chart(&t, &temps, true),
    ));

    // ── Boundary ────────────────────────────────────────────────────────────
    if samples.iter().any(|s| s.atoms_in.iter().chain(s.atoms_out.iter()).any(|&v| v != 0)) {
        let lines: Vec<(String, Vec<f64>)> = vec![
            ("atoms in".into(), samples.iter().map(|s| tot(&s.atoms_in)).collect()),
            ("atoms out".into(), samples.iter().map(|s| tot(&s.atoms_out)).collect()),
        ];
        o.push_str(&panel(
            "Matter across the boundary",
            "gross, never net. Equal and rising is a steady state; equal and flat is \
             a sealed box; diverging is accumulation.",
            &chart(&t, &lines, false),
        ));
    }

    o.push_str(
        "<p class=\"note\">These colours are arbitrary and mean nothing — they exist so \
         two lines can be told apart. This draws measurements of the world, not the world. \
         Rendering proper, when it comes, has to derive appearance from what things are \
         made of; nothing here does, and nothing here should be mistaken for it.</p>",
    );
    o
}

fn tot(v: &[i128]) -> f64 {
    v.iter().map(|&x| x as f64).sum()
}

fn panel(title: &str, note: &str, svg: &str) -> String {
    format!(
        "<div class=\"panel\"><h2>{}</h2><p>{}</p>{}</div>",
        esc(title),
        esc(note),
        svg
    )
}

/// Signed log, so a series spanning zero to 10¹⁸ stays legible and a negative
/// value (which should never happen, and therefore matters most) stays visible
/// instead of silently clamping to the floor.
fn slog(v: f64) -> f64 {
    if !v.is_finite() {
        return 0.0;
    }
    v.signum() * (1.0 + v.abs()).ln()
}

fn chart(t: &[f64], lines: &[(String, Vec<f64>)], log: bool) -> String {
    let f = |v: f64| if log { slog(v) } else { v };
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for (_, ys) in lines {
        for &y in ys {
            let y = f(y);
            if y.is_finite() {
                lo = lo.min(y);
                hi = hi.max(y);
            }
        }
    }
    if !lo.is_finite() || !hi.is_finite() {
        lo = 0.0;
        hi = 1.0;
    }
    if (hi - lo).abs() < 1e-12 {
        hi = lo + 1.0;
    }
    let (t0, t1) = (t.first().copied().unwrap_or(0.0), t.last().copied().unwrap_or(1.0));
    let tspan = if (t1 - t0).abs() < 1e-12 { 1.0 } else { t1 - t0 };

    let x = |v: f64| PAD_L + (v - t0) / tspan * (W - PAD_L - PAD_R);
    let y = |v: f64| PAD_T + (1.0 - (f(v) - lo) / (hi - lo)) * (H - PAD_T - PAD_B);

    let mut s = format!(r#"<svg viewBox="0 0 {} {}" width="100%">"#, W, H);
    // axes
    s.push_str(&format!(
        r#"<line class="ax" x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}"/>"#,
        PAD_L,
        H - PAD_B,
        W - PAD_R,
        H - PAD_B
    ));
    s.push_str(&format!(
        r#"<line class="ax" x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}"/>"#,
        PAD_L,
        PAD_T,
        PAD_L,
        H - PAD_B
    ));
    // y labels: the real values at the extremes, not the transformed ones
    let raw_hi = lines.iter().flat_map(|(_, v)| v.iter()).cloned().fold(f64::NEG_INFINITY, f64::max);
    let raw_lo = lines.iter().flat_map(|(_, v)| v.iter()).cloned().fold(f64::INFINITY, f64::min);
    s.push_str(&format!(
        r#"<text x="{:.1}" y="{:.1}" text-anchor="end">{}</text>"#,
        PAD_L - 8.0,
        PAD_T + 4.0,
        esc(&num(raw_hi))
    ));
    s.push_str(&format!(
        r#"<text x="{:.1}" y="{:.1}" text-anchor="end">{}</text>"#,
        PAD_L - 8.0,
        H - PAD_B,
        esc(&num(raw_lo))
    ));
    s.push_str(&format!(
        r#"<text x="{:.1}" y="{:.1}">tick {}</text>"#,
        PAD_L,
        H - PAD_B + 20.0,
        t0 as i64
    ));
    s.push_str(&format!(
        r#"<text x="{:.1}" y="{:.1}" text-anchor="end">tick {}</text>"#,
        W - PAD_R,
        H - PAD_B + 20.0,
        t1 as i64
    ));

    for (i, (name, ys)) in lines.iter().enumerate() {
        let c = SERIES[i % SERIES.len()];
        let mut d = String::new();
        for (k, &v) in ys.iter().enumerate() {
            let px = x(t.get(k).copied().unwrap_or(0.0));
            let py = y(v);
            if !py.is_finite() {
                continue;
            }
            d.push_str(&format!("{}{:.2} {:.2}", if d.is_empty() { "M" } else { "L" }, px, py));
            d.push(' ');
        }
        s.push_str(&format!(
            r#"<path d="{}" fill="none" stroke="{}" stroke-width="1.6" stroke-linejoin="round"/>"#,
            d, c
        ));
        let ly = PAD_T + 12.0 + i as f64 * 15.0;
        s.push_str(&format!(
            r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2.4"/>"#,
            W - PAD_R + 12.0,
            ly - 4.0,
            W - PAD_R + 28.0,
            ly - 4.0,
            c
        ));
        s.push_str(&format!(
            r#"<text x="{:.1}" y="{:.1}">{}</text>"#,
            W - PAD_R + 34.0,
            ly,
            esc(name)
        ));
    }
    s.push_str("</svg>");
    s
}

fn num(v: f64) -> String {
    if !v.is_finite() {
        return "—".into();
    }
    let a = v.abs();
    if a >= 1.0e5 || (a > 0.0 && a < 1.0e-2) {
        format!("{:.2e}", v)
    } else if a >= 1.0 {
        format!("{:.0}", v)
    } else {
        format!("{:.3}", v)
    }
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}
