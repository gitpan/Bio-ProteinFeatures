# $Id: test.pl,v 1.3 2003/09/07 00:22:00 cvspub Exp $
use Test::More qw(no_plan);
BEGIN{ use_ok( 'Bio::ProteinFeatures' ) };
ok(1);

$pf = new Bio::ProteinFeatures;

$sequence_string = 'TFDWRFAFTVEQGHREMIPVLPATMHGWIDQWVHSQACSRNNGDENCICPSLLM';

$pf->sequence($sequence_string);
ok($pf->sequence(), $sequence_string);

$f = $pf->features();
ok($f->{aa_prob}->{F} => 0.0555555555555555556);
ok($f->{first_order_entropy} => 2.89037175789616469);
ok($f->{transition_probability}->{neutral_polar} => 0.132075471698113208);

