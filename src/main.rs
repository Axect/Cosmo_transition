use peroxide::fuga::*;
use rayon::prelude::*;
use indicatif::{ProgressBar, ParallelProgressIterator};
use std::f64::consts::PI;

// Higgs mass
const MH: f64 = 125.0;
const LAMBDAC: f64 = 0.13;
const G: f64 = 0.65;
const GP: f64 = 0.35;
const YT: f64 = 1.0;
const NW: usize = 6;
const NZ: usize = 3;
#[allow(non_upper_case_globals)]
const Nt: usize = 12;
const CW: f64 = 5f64 / 6f64;
const CZ: f64 = 5f64 / 6f64;
#[allow(non_upper_case_globals)]
const Ct: f64 = 3f64 / 2f64;

#[allow(non_snake_case)]
fn main() {
    let phi_c_vec = linspace(0, 500, 10000);
    let V0_vec = phi_c_vec.fmap(V0);
    let V1loop_vec = phi_c_vec.fmap(V1loop);
    let V_eff_vec = V0_vec.add_v(&V1loop_vec);

    let mut plt = Plot2D::new();
    plt
        .set_domain(phi_c_vec)
        .insert_image(V0_vec)
        .insert_image(V1loop_vec)
        .insert_image(V_eff_vec.clone())
        .set_ylim((-5e8, 2e8))
        .set_style(PlotStyle::Nature)
        .tight_layout()
        .set_dpi(600)
        .set_path("one_loop.png")
        .savefig().unwrap();

    let T_vec = linspace(0.1, 200, 100);

    let result = T_vec.par_iter()
        .progress_with(ProgressBar::new(T_vec.len() as u64))
        .map(|&T| {
            let J_B_W = J_B(246f64, T, Boson::W);
            let J_B_Z = J_B(246f64, T, Boson::Z);
            let J_F_Top = J_F(246f64, T, Fermion::Top);
            ((J_B_W, J_B_Z), J_F_Top)
        })
        .collect::<Vec<_>>();

    let ((J_B_W, J_B_Z), J_F_Top): ((Vec<f64>, Vec<f64>), Vec<f64>) = result.into_iter().unzip();

    let V1loop_thermal_vec = T_vec.par_iter()
        .map(|&T| V1loop_thermal(246f64, T))
        .collect::<Vec<_>>();

    let mut plt = Plot2D::new();
    plt
        .set_domain(T_vec.clone())
        .insert_image(J_B_W)
        .insert_image(J_B_Z)
        .insert_image(J_F_Top)
        .set_style(PlotStyle::Nature)
        .set_legend(vec!["W", "Z", "Top"])
        .set_line_style(vec![(0, LineStyle::Dashed), (1, LineStyle::Dotted), (2, LineStyle::DashDot)])
        .tight_layout()
        .set_dpi(600)
        .set_path("J.png")
        .savefig().unwrap();

    let mut plt = Plot2D::new();
    plt
        .set_domain(T_vec)
        .insert_image(V1loop_thermal_vec)
        .set_style(PlotStyle::Nature)
        .set_dpi(600)
        .tight_layout()
        .set_path("one_loop_thermal_T.png")
        .savefig().unwrap();

    let phi_c_vec = linspace(0, 500, 100);
    let V0 = phi_c_vec.fmap(V0);
    let V1loop = phi_c_vec.fmap(V1loop);
    let V_eff = V0.add_v(&V1loop);
    let T = 100f64;
    let V1loop_thermal_vec = phi_c_vec.par_iter()
        .progress_with(ProgressBar::new(phi_c_vec.len() as u64))
        .map(|&phi_c| V1loop_thermal(phi_c, T))
        .collect::<Vec<_>>();
    let V1loop_thermal_vec = V1loop_thermal_vec.sub_s(V1loop_thermal_vec[0]);

    let V_tot = V_eff.add_v(&V1loop_thermal_vec);

    let one_loop_thermal_legend = format!(r"One-loop($T={}$)", T);

    let mut plt = Plot2D::new();
    plt
        .set_domain(phi_c_vec)
        .insert_image(V0.clone())
        .insert_image(V1loop)
        .insert_image(V1loop_thermal_vec)
        .insert_image(V_tot)
        .set_legend(vec!["Tree", "One-loop", one_loop_thermal_legend.as_str(), "Total"])
        .set_line_style(vec![(0, LineStyle::Dashed), (1, LineStyle::Dotted), (2, LineStyle::DashDot), (3, LineStyle::Solid)])
        .set_color(vec![(0, "red"), (1, "darkblue"), (2, "darkgreen"), (3, "black")])
        .set_style(PlotStyle::Nature)
        .set_dpi(600)
        .tight_layout()
        .set_path("one_loop_thermal_phi_c.png")
        .savefig().unwrap();

    let phi_c_vec = linspace(0, 500, 100);
    let T_vec = vec![50f64, 200f64, 300f64];
    let V1loop_thermal_vec = T_vec.par_iter()
        .progress_with(ProgressBar::new(T_vec.len() as u64))
        .map(|&T| phi_c_vec.clone().fmap(|phi_c| V1loop_thermal(phi_c, T)))
        .collect::<Vec<_>>();
    let V1loop_thermal_vec = V1loop_thermal_vec
        .into_iter()
        .map(|v| v.sub_s(v[0]))
        .collect::<Vec<_>>();
    let V_tot_vec = V1loop_thermal_vec
        .into_iter()
        .map(|v| v.add_v(&V_eff))
        .collect::<Vec<_>>();

    let mut plt = Plot2D::new();
    plt.set_domain(phi_c_vec);
    for V_tot in V_tot_vec {
        plt.insert_image(V_tot);
    }
    plt
        .set_legend(vec![r"$T = 50$", r"$T = 200$", r"$T = 300$"])
        .set_line_style(vec![(0, LineStyle::Dashed), (1, LineStyle::Dotted), (2, LineStyle::DashDot)])
        .set_color(vec![(0, "red"), (1, "darkblue"), (2, "darkgreen")])
        .set_style(PlotStyle::Nature)
        .set_xlabel(r"$\phi_c$")
        .set_ylabel(r"$V_{\rm tot}$")
        .set_dpi(600)
        .tight_layout()
        .set_path("V_tot.png")
        .savefig().unwrap();
}

#[allow(non_snake_case)]
pub fn m_W(phi_c: f64) -> f64 {
    (0.25 * G.powi(2) * phi_c.powi(2)).sqrt()
}

#[allow(non_snake_case)]
pub fn m_Z(phi_c: f64) -> f64 {
    (0.25 * (G.powi(2) + GP.powi(2)) * phi_c.powi(2)).sqrt()
}

pub fn m_t(phi_c: f64) -> f64 {
    (0.5 * YT.powi(2) * phi_c.powi(2)).sqrt()
}

#[allow(non_snake_case)]
pub fn V0(phi_c: f64) -> f64 {
    -0.5 * MH.powi(2) * phi_c.powi(2) + 0.25 * LAMBDAC * phi_c.powi(4)
}

#[allow(non_snake_case)]
pub fn V1loop(phi_c: f64) -> f64 {
    let m_W = m_W(phi_c);
    let m_Z = m_Z(phi_c);
    let m_t = m_t(phi_c);
    1f64 / (64f64 * PI.powi(2)) * (
        NW as f64 * m_W.powi(4) * ((m_W.powi(2) / MH.powi(2)).ln() - CW)
        + NZ as f64 * m_Z.powi(4) * ((m_Z.powi(2) / MH.powi(2)).ln() - CZ)
        - Nt as f64 * m_t.powi(4) * ((m_t.powi(2) / MH.powi(2)).ln() - Ct)
    )
}

#[derive(Debug, Copy, Clone)]
pub enum Boson {
    W,
    Z
}

#[derive(Debug, Copy, Clone)]
pub enum Fermion {
    Top,
}

#[allow(non_snake_case)]
pub fn J_B(phi_c: f64, T: f64, boson: Boson) -> f64 {
    let m_B_T = match boson {
        Boson::W => m_W(phi_c).powi(2) / T.powi(2),
        Boson::Z => m_Z(phi_c).powi(2) / T.powi(2),
    };
    let lambda = 100f64;
    let f = |t: f64| {
        t.powi(2) * (1f64 - (-((lambda * t).powi(2) + m_B_T).sqrt()).exp()).ln()
    };
    lambda.powi(3) * integrate(f, (0f64, 1f64), G7K15R(1e-3, 15))
}

#[allow(non_snake_case)]
pub fn J_F(phi_c: f64, T: f64, fermion: Fermion) -> f64 {
    let m_F = match fermion {
        Fermion::Top => m_t(phi_c).powi(2) / T.powi(2),
    };
    let lambda = 100f64;
    let f = |t: f64| {
        t.powi(2) * (1f64 + (-((lambda * t).powi(2) + m_F).sqrt()).exp()).ln()
    };
    lambda.powi(3) * integrate(f, (0f64, 1f64), G7K15R(1e-3, 15))
}

#[allow(non_snake_case)]
pub fn V1loop_thermal(phi_c: f64, T: f64) -> f64 {
    let J_W = J_B(phi_c, T, Boson::W);
    let J_Z = J_B(phi_c, T, Boson::Z);
    let J_Top = J_F(phi_c, T, Fermion::Top);
    T.powi(4) / (2f64 * PI.powi(2)) * (NW as f64 * J_W + NZ as f64 * J_Z - Nt as f64 * J_Top)
}