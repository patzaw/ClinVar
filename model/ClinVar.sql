-- MySQL Script generated by MySQL Workbench
-- Wed Apr  4 11:40:12 2018
-- Model: New Model    Version: 1.0
-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL,ALLOW_INVALID_DATES';

-- -----------------------------------------------------
-- Schema ClinVar
-- -----------------------------------------------------

-- -----------------------------------------------------
-- Schema ClinVar
-- -----------------------------------------------------
CREATE SCHEMA IF NOT EXISTS `ClinVar` DEFAULT CHARACTER SET utf8 ;
USE `ClinVar` ;

-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_clinSigOrder`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_clinSigOrder` (
  `label` VARCHAR(45) NOT NULL,
  `order` INT NOT NULL COMMENT 'An integer indicating the order of clinical significance: the lowest the best clinical outcome.',
  PRIMARY KEY (`label`))
ENGINE = InnoDB
COMMENT = 'Not provided by ClinVar. Manually generated.';


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_revStatOrder`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_revStatOrder` (
  `label` VARCHAR(45) NOT NULL,
  `order` INT NOT NULL COMMENT 'An integer indicating the order of review status: the lowest hte least significant.',
  PRIMARY KEY (`label`))
ENGINE = InnoDB
COMMENT = 'Not provided by ClinVar. Manually generated.';


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_ReferenceClinVarAssertion`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_ReferenceClinVarAssertion` (
  `cvs` INT NOT NULL,
  `id` INT NOT NULL,
  `accession` VARCHAR(45) NOT NULL,
  `assertion` VARCHAR(45) NULL,
  `reviewStatus` VARCHAR(45) NOT NULL,
  `clinicalSignificance` VARCHAR(45) NOT NULL,
  `explanation` VARCHAR(45) NULL,
  `title` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE INDEX `cvs_UNIQUE` (`cvs` ASC),
  INDEX `fk_ClinVar_ReferenceClinVarAssertion_ClinVar_clinSigOrder1_idx` (`clinicalSignificance` ASC),
  CONSTRAINT `fk_ClinVar_ReferenceClinVarAssertion_ClinVar_clinSigOrder1`
    FOREIGN KEY (`clinicalSignificance`)
    REFERENCES `ClinVar`.`ClinVar_clinSigOrder` (`label`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_ClinVar_ReferenceClinVarAssertion_ClinVar_revStatOrder1`
    FOREIGN KEY (`reviewStatus`)
    REFERENCES `ClinVar`.`ClinVar_revStatOrder` (`label`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_rcvaInhMode`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_rcvaInhMode` (
  `rcvaId` INT NOT NULL,
  `inhMode` VARCHAR(45) NOT NULL,
  CONSTRAINT `fk_cv_rcvaSubmitters_cv_ReferenceClinVarAssertion0`
    FOREIGN KEY (`rcvaId`)
    REFERENCES `ClinVar`.`ClinVar_ReferenceClinVarAssertion` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB
COMMENT = 'A few RCVA have several mode of inheritance. These cases are handled by this table.';


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_rcvaObservedIn`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_rcvaObservedIn` (
  `rcvaId` INT NOT NULL,
  `origin` VARCHAR(45) NOT NULL,
  `taxonomyId` VARCHAR(45) NOT NULL,
  `species` VARCHAR(45) NOT NULL,
  `affectedStatus` VARCHAR(45) NOT NULL,
  `numberTested` VARCHAR(45) NULL,
  `sampleDescription` VARCHAR(45) NULL,
  `nbObsInTraitSet` INT NOT NULL,
  `nbSampleTraitSet` INT NOT NULL,
  CONSTRAINT `fk_cv_rcvaSubmitters_cv_ReferenceClinVarAssertion00`
    FOREIGN KEY (`rcvaId`)
    REFERENCES `ClinVar`.`ClinVar_ReferenceClinVarAssertion` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_ClinVarAssertions`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_ClinVarAssertions` (
  `cvs` INT NOT NULL,
  `id` INT NOT NULL,
  `accession` VARCHAR(45) NOT NULL,
  `clinicalSignificance` VARCHAR(45) NULL,
  PRIMARY KEY (`id`),
  INDEX `fk_ClinVar_ClinVarAssertions_ClinVar_ReferenceClinVarAssert_idx` (`cvs` ASC),
  CONSTRAINT `fk_ClinVar_ClinVarAssertions_ClinVar_ReferenceClinVarAssertion1`
    FOREIGN KEY (`cvs`)
    REFERENCES `ClinVar`.`ClinVar_ReferenceClinVarAssertion` (`cvs`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_cvaObservedIn`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_cvaObservedIn` (
  `cvaId` INT NOT NULL,
  `origin` VARCHAR(45) NOT NULL,
  `species` VARCHAR(45) NOT NULL,
  `affectedStatus` VARCHAR(45) NOT NULL,
  CONSTRAINT `fk_ClinVar_cvaObservedIn_ClinVar_ClinVarAssertions1`
    FOREIGN KEY (`cvaId`)
    REFERENCES `ClinVar`.`ClinVar_ClinVarAssertions` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_cvaSubmitters`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_cvaSubmitters` (
  `cvaId` INT NOT NULL,
  `submitter` VARCHAR(45) NOT NULL,
  `primary` TINYINT NOT NULL,
  CONSTRAINT `fk_ClinVar_cvaObservedIn_ClinVar_ClinVarAssertions10`
    FOREIGN KEY (`cvaId`)
    REFERENCES `ClinVar`.`ClinVar_ClinVarAssertions` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_sourceFiles`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_sourceFiles` (
  `url` VARCHAR(45) NOT NULL,
  `current` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`url`))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_variants`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_variants` (
  `id` INT NOT NULL,
  `type` VARCHAR(45) NOT NULL,
  `name` VARCHAR(45) NULL,
  PRIMARY KEY (`id`))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_rcvaVariant`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_rcvaVariant` (
  `varId` INT NOT NULL,
  `rcvaId` INT NOT NULL,
  INDEX `fk_ClinVar_rcvaVariant_ClinVar_ReferenceClinVarAssertion1_idx` (`rcvaId` ASC),
  CONSTRAINT `fk_ClinVar_rcvaVariant_ClinVar_ReferenceClinVarAssertion1`
    FOREIGN KEY (`rcvaId`)
    REFERENCES `ClinVar`.`ClinVar_ReferenceClinVarAssertion` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_ClinVar_rcvaVariant_ClinVar_variants1`
    FOREIGN KEY (`varId`)
    REFERENCES `ClinVar`.`ClinVar_variants` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_varCytoLoc`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_varCytoLoc` (
  `varId` INT NOT NULL,
  `location` VARCHAR(45) NOT NULL,
  INDEX `fk_ClinVar_varCytoLoc_ClinVar_variants1_idx` (`varId` ASC),
  CONSTRAINT `fk_ClinVar_varCytoLoc_ClinVar_variants1`
    FOREIGN KEY (`varId`)
    REFERENCES `ClinVar`.`ClinVar_variants` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_varAttributes`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_varAttributes` (
  `varId` INT NOT NULL,
  `Type` VARCHAR(45) NOT NULL,
  `integerValue` INT NULL,
  `Change` VARCHAR(45) NULL,
  `value` VARCHAR(45) NOT NULL,
  INDEX `fk_ClinVar_varCytoLoc_ClinVar_variants1_idx` (`varId` ASC),
  CONSTRAINT `fk_ClinVar_varCytoLoc_ClinVar_variants10`
    FOREIGN KEY (`varId`)
    REFERENCES `ClinVar`.`ClinVar_variants` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_varXRef`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_varXRef` (
  `varId` INT NOT NULL,
  `id` VARCHAR(45) NOT NULL,
  `db` VARCHAR(45) NOT NULL,
  `type` VARCHAR(45) NULL,
  INDEX `fk_ClinVar_varCytoLoc_ClinVar_variants1_idx` (`varId` ASC),
  CONSTRAINT `fk_ClinVar_varCytoLoc_ClinVar_variants11`
    FOREIGN KEY (`varId`)
    REFERENCES `ClinVar`.`ClinVar_variants` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_varNames`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_varNames` (
  `varId` INT NOT NULL,
  `name` VARCHAR(45) NOT NULL,
  `type` VARCHAR(45) NOT NULL,
  INDEX `fk_ClinVar_varCytoLoc_ClinVar_variants1_idx` (`varId` ASC),
  CONSTRAINT `fk_ClinVar_varCytoLoc_ClinVar_variants12`
    FOREIGN KEY (`varId`)
    REFERENCES `ClinVar`.`ClinVar_variants` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_entrezNames`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_entrezNames` (
  `entrez` INT NOT NULL,
  `name` VARCHAR(45) NOT NULL,
  `symbol` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`entrez`))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_varEntrez`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_varEntrez` (
  `varId` INT NOT NULL,
  `entrez` INT NOT NULL,
  `type` VARCHAR(45) NOT NULL,
  INDEX `fk_ClinVar_varCytoLoc_ClinVar_variants1_idx` (`varId` ASC),
  CONSTRAINT `fk_ClinVar_varCytoLoc_ClinVar_variants13`
    FOREIGN KEY (`varId`)
    REFERENCES `ClinVar`.`ClinVar_variants` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_ClinVar_varEntrez_ClinVar_entrezNames1`
    FOREIGN KEY (`varId`)
    REFERENCES `ClinVar`.`ClinVar_entrezNames` (`entrez`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_varSeqLoc`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_varSeqLoc` (
  `varId` INT NOT NULL,
  `Accession` VARCHAR(45) NOT NULL,
  `alternateAllele` VARCHAR(45) NULL,
  `Assembly` VARCHAR(45) NOT NULL,
  `AssemblyAccessionVersion` VARCHAR(45) NOT NULL,
  `AssemblyStatus` VARCHAR(45) NOT NULL,
  `Chr` VARCHAR(45) NOT NULL,
  `display_start` INT NOT NULL,
  `display_stop` INT NOT NULL,
  `innerStart` INT NULL,
  `innerStop` INT NULL,
  `outerStart` INT NULL,
  `outerStop` INT NULL,
  `referenceAllele` VARCHAR(45) NULL,
  `start` INT NULL,
  `stop` INT NULL,
  `Strand` VARCHAR(45) NULL,
  `variantLength` INT NULL,
  INDEX `fk_ClinVar_varCytoLoc_ClinVar_variants1_idx` (`varId` ASC),
  CONSTRAINT `fk_ClinVar_varCytoLoc_ClinVar_variants14`
    FOREIGN KEY (`varId`)
    REFERENCES `ClinVar`.`ClinVar_variants` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_traits`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_traits` (
  `id` INT NOT NULL,
  `name` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = big5
COMMENT = '		';


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_rcvaTraits`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_rcvaTraits` (
  `rcvaId` INT NOT NULL,
  `t.id` INT NOT NULL,
  `traitType` VARCHAR(45) NOT NULL,
  INDEX `fk_ClinVar_rcvaTraits_ClinVar_traits1_idx` (`t.id` ASC),
  CONSTRAINT `fk_ClinVar_rcvaTraits_ClinVar_traits1`
    FOREIGN KEY (`t.id`)
    REFERENCES `ClinVar`.`ClinVar_traits` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_ClinVar_rcvaTraits_ClinVar_ReferenceClinVarAssertion1`
    FOREIGN KEY (`rcvaId`)
    REFERENCES `ClinVar`.`ClinVar_ReferenceClinVarAssertion` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_traitCref`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_traitCref` (
  `t.id` INT NOT NULL,
  `id` VARCHAR(45) NOT NULL,
  `db` VARCHAR(45) NOT NULL,
  `type` VARCHAR(45) NULL,
  CONSTRAINT `fk_ClinVar_traitCref_ClinVar_traits1`
    FOREIGN KEY (`t.id`)
    REFERENCES `ClinVar`.`ClinVar_traits` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `ClinVar`.`ClinVar_traitNames`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `ClinVar`.`ClinVar_traitNames` (
  `t.id` INT NOT NULL,
  `name` VARCHAR(45) NOT NULL,
  `type` VARCHAR(45) NOT NULL COMMENT '\"Alternate\" or \"Preferred\"',
  CONSTRAINT `fk_ClinVar_traitCref_ClinVar_traits10`
    FOREIGN KEY (`t.id`)
    REFERENCES `ClinVar`.`ClinVar_traits` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;
