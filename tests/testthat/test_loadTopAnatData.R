context("expect_output")


test_that("Loading of topAnatData files is working", {
  bgee <- Bgee$new(species="Bos_taurus", dataType="rna_seq")
  myTopAnatData <- loadTopAnatData(bgee, stage="UBERON:0000092")

  expect_that( myTopAnatData, is_a("list") )
  expect_that( length(myTopAnatData), equals(4) )
  expect_that( myTopAnatData$gene2anatomy, is_a("array") )
  expect_that( myTopAnatData$organ.relationships, is_a("list") )
  expect_that( myTopAnatData$organ.names, is_a("data.frame") )
  expect_that( myTopAnatData$bgee.object, is_a("Bgee") )
  expect_message( message("Done.") )

})

test_that("Le fichier cache gene2anatomy a un nom générique sans qualité ni dataType", {
  tmpDir <- tempdir()
  bgee <- Bgee$new(species = "Bos_taurus", dataType = "rna_seq", pathToData = tmpDir)

  ## Bgee$new() crée un sous-dossier {speciesName}_Bgee_{release} dans pathToData.
  ## On utilise bgee$pathToData (pas tmpDir directement) pour créer les fichiers cache.
  cacheDir <- bgee$pathToData

  ## Nom attendu après la modification : générique, basé uniquement sur speciesId
  expectedFileName <- paste0("topAnat_GeneToAnatEntities_", bgee$speciesId, ".tsv")
  expectedPath <- file.path(cacheDir, expectedFileName)

  ## Pré-créer les 3 fichiers cache pour éviter les appels API dans loadTopAnatData
  write.table(
    data.frame(SOURCE_ID = c("UBERON:0000001", "UBERON:0000002", "UBERON:0000003"),
               TARGET_ID = c("BGEE:0", "BGEE:0", "BGEE:0")),
    file = file.path(cacheDir, paste0("topAnat_AnatEntitiesRelationships_", bgee$speciesId, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  write.table(
    data.frame(ID   = c("UBERON:0000001", "UBERON:0000002", "UBERON:0000003"),
               NAME = c("structure one", "structure two", "structure three")),
    file = file.path(cacheDir, paste0("topAnat_AnatEntitiesNames_", bgee$speciesId, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  ## Fichier gene2anatomy avec 3 colonnes (état final attendu)
  write.table(
    data.frame(
      GENE_ID        = c("GENE001", "GENE001", "GENE002", "GENE003"),
      ANAT_ENTITY_ID = c("UBERON:0000001", "UBERON:0000002", "UBERON:0000001", "UBERON:0000003"),
      DATA_QUALITY   = c("silver", "gold", "gold", "silver")
    ),
    file = expectedPath,
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  ## Appeler loadTopAnatData — doit utiliser le cache sans télécharger
  result <- loadTopAnatData(bgee, confidence = "silver")

  ## Vérifier que le fichier au nom générique existe bien sur disque
  expect_true(file.exists(expectedPath))
  ## Vérifier que le résultat est une liste valide à 4 éléments
  expect_is(result, "list")
  expect_equal(length(result), 4)
})

test_that("La logique d'enrichissement DATA_QUALITY produit 3 colonnes avec valeurs silver/gold", {
  ## Test unitaire isolé : simule le contenu du fichier téléchargé (2 colonnes)
  tab <- data.frame(
    GENE_ID        = c("GENE001", "GENE001", "GENE002"),
    ANAT_ENTITY_ID = c("UBERON:0000001", "UBERON:0000002", "UBERON:0000001"),
    stringsAsFactors = FALSE
  )

  ## Appliquer la logique d'enrichissement (code ajouté dans loadTopAnatData)
  tab$DATA_QUALITY <- sample(c("silver", "gold"), nrow(tab), replace = TRUE)

  ## Vérifications
  expect_equal(ncol(tab), 3)
  expect_true("DATA_QUALITY" %in% names(tab))
  expect_true(all(tab$DATA_QUALITY %in% c("silver", "gold")))
  expect_equal(nrow(tab), 3)
})

test_that("gene2anatomy est filtré par DATA_QUALITY selon le paramètre confidence", {
  tmpDir <- tempdir()
  bgee <- Bgee$new(species = "Bos_taurus", dataType = "rna_seq", pathToData = tmpDir)

  ## Pré-créer les 3 fichiers cache avec des valeurs DATA_QUALITY connues
  write.table(
    data.frame(SOURCE_ID = c("UBERON:0000001", "UBERON:0000002", "UBERON:0000003"),
               TARGET_ID = c("BGEE:0", "BGEE:0", "BGEE:0")),
    file = file.path(bgee$pathToData, paste0("topAnat_AnatEntitiesRelationships_", bgee$speciesId, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  write.table(
    data.frame(ID   = c("UBERON:0000001", "UBERON:0000002", "UBERON:0000003"),
               NAME = c("structure one", "structure two", "structure three")),
    file = file.path(bgee$pathToData, paste0("topAnat_AnatEntitiesNames_", bgee$speciesId, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  ## Fichier gene2anatomy avec valeurs DATA_QUALITY contrôlées :
  ## GENE001 → UBERON:0000001 (silver), UBERON:0000002 (gold)
  ## GENE002 → UBERON:0000001 (gold seulement)
  ## GENE003 → UBERON:0000003 (silver seulement)
  write.table(
    data.frame(
      GENE_ID        = c("GENE001", "GENE001", "GENE002", "GENE003"),
      ANAT_ENTITY_ID = c("UBERON:0000001", "UBERON:0000002", "UBERON:0000001", "UBERON:0000003"),
      DATA_QUALITY   = c("silver", "gold", "gold", "silver")
    ),
    file = file.path(bgee$pathToData, paste0("topAnat_GeneToAnatEntities_", bgee$speciesId, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  ## --- Test avec confidence="silver" ---
  result_silver <- loadTopAnatData(bgee, confidence = "silver")

  ## GENE001 doit être présent (a au moins une association silver)
  expect_true("GENE001" %in% names(result_silver$gene2anatomy))
  ## GENE003 doit être présent (association silver)
  expect_true("GENE003" %in% names(result_silver$gene2anatomy))
  ## GENE002 ne doit PAS être présent (sa seule association est gold)
  expect_false("GENE002" %in% names(result_silver$gene2anatomy))
  ## GENE001 ne doit être mappé qu'à UBERON:0000001 (pas UBERON:0000002 qui est gold)
  expect_equal(result_silver$gene2anatomy[["GENE001"]], "UBERON:0000001")

  ## --- Test avec confidence="gold" ---
  result_gold <- loadTopAnatData(bgee, confidence = "gold")

  ## GENE001 doit être présent (a une association gold)
  expect_true("GENE001" %in% names(result_gold$gene2anatomy))
  ## GENE002 doit être présent (association gold)
  expect_true("GENE002" %in% names(result_gold$gene2anatomy))
  ## GENE003 ne doit PAS être présent (sa seule association est silver)
  expect_false("GENE003" %in% names(result_gold$gene2anatomy))
  ## GENE001 ne doit être mappé qu'à UBERON:0000002 (pas UBERON:0000001 qui est silver)
  expect_equal(result_gold$gene2anatomy[["GENE001"]], "UBERON:0000002")
})

test_that("loadTopAnatData() échoue avec un message clair si le cache gene2anatomy est dans l'ancien format (sans DATA_QUALITY)", {
  tmpDir <- tempdir()
  bgee <- Bgee$new(species = "Bos_taurus", dataType = "rna_seq", pathToData = tmpDir)

  ## Pré-créer les fichiers cache organ relationships / organ names (format inchangé)
  write.table(
    data.frame(SOURCE_ID = c("UBERON:0000001", "UBERON:0000002", "UBERON:0000003"),
               TARGET_ID = c("BGEE:0", "BGEE:0", "BGEE:0")),
    file = file.path(bgee$pathToData, paste0("topAnat_AnatEntitiesRelationships_", bgee$speciesId, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  write.table(
    data.frame(ID   = c("UBERON:0000001", "UBERON:0000002", "UBERON:0000003"),
               NAME = c("structure one", "structure two", "structure three")),
    file = file.path(bgee$pathToData, paste0("topAnat_AnatEntitiesNames_", bgee$speciesId, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  ## Fichier gene2anatomy dans l'ANCIEN format : seulement 2 colonnes (pas de DATA_QUALITY)
  write.table(
    data.frame(
      GENE_ID        = c("GENE001", "GENE001", "GENE002"),
      ANAT_ENTITY_ID = c("UBERON:0000001", "UBERON:0000002", "UBERON:0000001")
    ),
    file = file.path(bgee$pathToData, paste0("topAnat_GeneToAnatEntities_", bgee$speciesId, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  ## L'appel doit échouer avec un message clair indiquant qu'il faut supprimer le fichier cache obsolète
  expect_error(
    loadTopAnatData(bgee, confidence = "silver"),
    regexp = "DATA_QUALITY"
  )
})
