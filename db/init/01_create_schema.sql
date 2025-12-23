DROP TABLE IF EXISTS `local_drugs`;
CREATE TABLE `local_drugs`  (
  `id` int NOT NULL AUTO_INCREMENT,
  `normalized_name` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_uca1400_ai_ci NOT NULL,
  `local_name` varchar(500) CHARACTER SET utf8mb4 COLLATE utf8mb4_uca1400_ai_ci NULL DEFAULT NULL,
  `chembl_name` varchar(500) CHARACTER SET utf8mb4 COLLATE utf8mb4_uca1400_ai_ci NULL DEFAULT NULL,
  `molecule_chembl_id` varchar(50) CHARACTER SET utf8mb4 COLLATE utf8mb4_uca1400_ai_ci NULL DEFAULT NULL,
  `smiles` text CHARACTER SET utf8mb4 COLLATE utf8mb4_uca1400_ai_ci NULL,
  `conformer_source` varchar(10) CHARACTER SET utf8mb4 COLLATE utf8mb4_uca1400_ai_ci NULL DEFAULT NULL,
  `molblock` mediumtext CHARACTER SET utf8mb4 COLLATE utf8mb4_uca1400_ai_ci NULL,
  PRIMARY KEY (`id`) USING BTREE,
  UNIQUE INDEX `uq_iran_drugs_normalized_name`(`normalized_name` ASC) USING BTREE,
  INDEX `idx_iran_drugs_chembl_id`(`molecule_chembl_id` ASC) USING BTREE
) ENGINE = InnoDB AUTO_INCREMENT = 764 CHARACTER SET = utf8mb4 COLLATE = utf8mb4_uca1400_ai_ci ROW_FORMAT = Dynamic;

CREATE TABLE chembl_approved_drugs (
    id INT NOT NULL AUTO_INCREMENT,
    molecule_chembl_id VARCHAR(50) NOT NULL UNIQUE,
    name VARCHAR(255),
    smiles TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (id)
);
SET FOREIGN_KEY_CHECKS = 1;

CREATE TABLE chembl_sync_state (
    id INT PRIMARY KEY,
    last_offset INT NOT NULL DEFAULT 0,
    is_completed BOOLEAN NOT NULL DEFAULT FALSE,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP
);


CREATE TABLE local_drug_import_log (
    id BIGINT AUTO_INCREMENT PRIMARY KEY,

    normalized_name VARCHAR(255) NOT NULL,
    raw_name VARCHAR(255) NULL,

    chembl_id VARCHAR(50) NULL,
    chembl_name VARCHAR(255) NULL,
    smiles TEXT NULL,

    status ENUM(
        'ALREADY_EXISTS',
        'NOT_IN_CHEMBL',
        'INVALID_SMILES',
        'TOO_LARGE',
        'CONFORMER_FAILED',
        'EXCLUDED'
    ) NOT NULL,

    message TEXT NULL,

    num_atoms INT NULL,

    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);


CREATE TABLE users (
    id BIGINT AUTO_INCREMENT PRIMARY KEY,
    username VARCHAR(100) NOT NULL UNIQUE,
    password_hash VARCHAR(255) NOT NULL,
    is_active BOOLEAN DEFAULT TRUE,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

INSERT INTO users (username, password_hash)
VALUES (
  'admin',
  SHA2('admin123', 256)
);

CREATE TABLE chembl_targets (
    target_chembl_id   VARCHAR(32) NOT NULL,
    target_name        VARCHAR(255),
    organism           VARCHAR(128),
    target_type        VARCHAR(64),

    PRIMARY KEY (target_chembl_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;


CREATE TABLE proteins (
    uniprot_id     VARCHAR(16) NOT NULL,
    protein_name   VARCHAR(255),
    gene_name      VARCHAR(64),
    organism       VARCHAR(128),

    PRIMARY KEY (uniprot_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;


CREATE TABLE drug_targets (
    molecule_chembl_id  VARCHAR(32) NOT NULL,
    target_chembl_id    VARCHAR(32) NOT NULL,

    PRIMARY KEY (molecule_chembl_id, target_chembl_id),

    CONSTRAINT fk_dt_drug
        FOREIGN KEY (molecule_chembl_id)
        REFERENCES chembl_approved_drugs (molecule_chembl_id)
        ON DELETE CASCADE,

    CONSTRAINT fk_dt_target
        FOREIGN KEY (target_chembl_id)
        REFERENCES chembl_targets (target_chembl_id)
        ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;


CREATE TABLE target_proteins (
    target_chembl_id  VARCHAR(32) NOT NULL,
    uniprot_id        VARCHAR(16) NOT NULL,

    PRIMARY KEY (target_chembl_id, uniprot_id),

    CONSTRAINT fk_tp_target
        FOREIGN KEY (target_chembl_id)
        REFERENCES chembl_targets (target_chembl_id)
        ON DELETE CASCADE,

    CONSTRAINT fk_tp_protein
        FOREIGN KEY (uniprot_id)
        REFERENCES proteins (uniprot_id)
        ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;


CREATE TABLE reactome_pathways (
    pathway_id        VARCHAR(32) NOT NULL,
    pathway_name      VARCHAR(255) NOT NULL,
    top_level_class   VARCHAR(128),
    species           VARCHAR(64),

    PRIMARY KEY (pathway_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE protein_pathways (
    uniprot_id   VARCHAR(16) NOT NULL,
    pathway_id   VARCHAR(32) NOT NULL,

    PRIMARY KEY (uniprot_id, pathway_id),

    CONSTRAINT fk_pp_protein
        FOREIGN KEY (uniprot_id)
        REFERENCES proteins (uniprot_id)
        ON DELETE CASCADE,

    CONSTRAINT fk_pp_pathway
        FOREIGN KEY (pathway_id)
        REFERENCES reactome_pathways (pathway_id)
        ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;


CREATE TABLE drug_pathways (
    molecule_chembl_id  VARCHAR(32) NOT NULL,
    pathway_id          VARCHAR(32) NOT NULL,

    PRIMARY KEY (molecule_chembl_id, pathway_id),

    CONSTRAINT fk_dp_drug
        FOREIGN KEY (molecule_chembl_id)
        REFERENCES chembl_approved_drugs (molecule_chembl_id)
        ON DELETE CASCADE,

    CONSTRAINT fk_dp_pathway
        FOREIGN KEY (pathway_id)
        REFERENCES reactome_pathways (pathway_id)
        ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE protein_reactome_status (
    uniprot_id     VARCHAR(16) PRIMARY KEY,
    has_pathway    BOOLEAN NOT NULL,
    checked_at     DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP
);

ALTER TABLE local_drugs
ADD COLUMN targets_synced TINYINT DEFAULT 0;

CREATE INDEX idx_dt_target ON drug_targets (target_chembl_id);
CREATE INDEX idx_tp_protein ON target_proteins (uniprot_id);
CREATE INDEX idx_pp_pathway ON protein_pathways (pathway_id);
CREATE INDEX idx_dp_pathway ON drug_pathways (pathway_id);
