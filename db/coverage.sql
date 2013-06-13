CREATE TABLE coverage (
    sample_id  TEXT    NOT NULL,
    chromosome TEXT    NOT NULL,
    position   INTEGER NOT NULL,
    gap_cov    INTEGER NOT NULL,
    nogap_cov  INTEGER NOT NULL,
    PRIMARY KEY ( sample_id, chromosome, position )
);
