// Reference counting helper
// TODO: delete me and replace with new/drop or `Rc`
#[macro_export]
macro_rules! refcnt {
    ($count:expr, $delta:expr) => {{
        let rc = $count + $delta;
        $count = rc;
        debug_assert!(rc >= 0);
        rc
    }};
}
