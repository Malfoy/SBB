#[allow(dead_code)]
mod utils {
    include!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/vendor/bloomybloom/src/utils.rs"
    ));
}

#[allow(dead_code)]
mod decyclers {
    include!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/vendor/bloomybloom/src/decyclers.rs"
    ));
}

#[allow(dead_code)]
mod minimizers {
    include!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/vendor/bloomybloom/src/minimizers.rs"
    ));
}

#[allow(dead_code)]
mod bloomybloom_vendor;
#[allow(dead_code)]
mod bloomybloom_adapter;

pub mod bloom;
pub mod classify;
pub mod fastx;
pub mod writer;
