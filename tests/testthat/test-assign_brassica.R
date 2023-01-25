test_that('brassica ASV assigned to multiple varieties', {
  brassica <- 'ATCCTGGGTTACGCGAACAAAACAGAGTTTAGAAAGCGG'
  ref <- test_path('..',
                   'testdata',
                   'trnLGH.fasta')
  expect_gt(object = nrow(assignSpecies_mod(brassica, ref)),
            expected = 1)
})
