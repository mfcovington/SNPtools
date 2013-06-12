CREATE TABLE coverage (
    chromosome  TEXT    NOT NULL,
    position    INTEGER NOT NULL,
    gap_depth   INTEGER NOT NULL,
    nogap_depth INTEGER NOT NULL,
    PRIMARY KEY (chromosome, position)
);
