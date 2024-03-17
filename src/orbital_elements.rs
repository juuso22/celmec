/// Keplerian elements of a 2-body system
///
/// **Fields**:
///
/// e = eccentricity
///
/// longitude_of_the_ascending_node = longitude of the ascending node
///
/// tau = perihelion time
///
/// a = semi-major axis
///
/// i = inclination
///
/// omega = argument of perihelion
pub struct KeplerianElements {
    pub e: f64,
    pub longitude_of_the_ascending_node: f64,
    pub tau: f64,
    pub a: f64,
    pub iota: f64,
    pub omega: f64,
}
