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
