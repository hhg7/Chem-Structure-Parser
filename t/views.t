#!/usr/bin/env perl
# The flat views over a parsed structure, and the summary.
require 5.010;
use strict;
use warnings FATAL => 'all';
use Cwd 'abs_path';
use File::Basename 'dirname';
use Scalar::Util 'refaddr';
use Chem::Structure::Parser;
use Test::Exception;
use Test::More;

my $data = dirname(abs_path(__FILE__)) . '/data';
my $i = structure_info("$data/mini.pdb");

#--------
# structure_atoms -- the shape to hand to a distance calculation
#--------
{
	my $atoms = structure_atoms($i);
	is(ref $atoms, 'ARRAY', 'structure_atoms: an array reference');
	is(scalar @$atoms, 65,
		'one entry per distinct atom: the two conformers of one CB are one atom');
	is($atoms->[0]{name}, 'N', 'in file order');
	is($atoms->[0]{chain}, 'A', 'each atom carries its chain');
	is($atoms->[0]{resname}, 'MET', 'and its residue name');
	is($atoms->[0]{resseq}, 1, 'and its residue number');
	is($atoms->[0]{reskey}, '1', 'and the key it is filed under');
	ok(defined $atoms->[0]{x}, 'and its coordinates');

	my $a = structure_atoms($i, 'B');
	is(scalar @$a, 12, 'structure_atoms: one chain on request');
	is($a->[0]{chain}, 'B', 'and it is that chain');

	# the copies must not be the residue-s own atom hashes: writing to one of
	# these should not scribble on the structure
	$a->[0]{x} = 'scribbled';
	isnt($i->{chains}{B}{residues}{1}{atoms}{P}{x}, 'scribbled',
		'structure_atoms returns copies, not the hashes inside the structure');

	throws_ok { structure_atoms($i, 'Z') } qr/no chain 'Z'/, 'a chain that is not there dies';
}

#--------
# structure_residues -- the same residues that are in the nested structure,
# not copies, so that walking them and looking one up agree
#--------
{
	my $res = structure_residues($i);
	is(scalar @$res, 17, 'structure_residues: every residue in the file');
	is($res->[0]{resname}, 'MET', 'in file order');
	is(refaddr($res->[0]), refaddr($i->{chains}{A}{residues}{1}),
		'and they are the very residues in the structure, not copies of them');
	is(scalar @{ structure_residues($i, 'B') }, 4, 'one chain on request');
}

#--------
# structure_ligands
#--------
{
	my $lig = structure_ligands($i);
	is_deeply([ sort keys %$lig ], [ 'NAG_A_201', 'ZN_A_202' ],
		'structure_ligands: the heterogens that are not water, keyed by name, chain and number');
	is($lig->{NAG_A_201}{resname}, 'NAG', 'and the residue itself is the value');
	ok(!grep({ /HOH/ } keys %$lig), 'water is not a ligand');
}

#--------
# sequences
#--------
{
	is_deeply(structure_sequences($i), { A => 'MAGCMHHSC', B => 'ACGT' }, 'structure_sequences');
	throws_ok { chain_sequence($i, 'Z') } qr/no chain 'Z'/, 'chain_sequence: unknown chain dies';
	throws_ok { chain_sequence($i) } qr/no chain given/, 'chain_sequence: no chain dies';
	throws_ok { chain_sequence($i, 'A', 'guessed') } qr/observed.*seqres/,
		'chain_sequence: an unknown kind of sequence dies';
}

#--------
# structure_summary
#--------
{
	my $s = structure_summary($i);
	like($s, qr/9XYZ/, 'summary: the id');
	like($s, qr/A SMALL TEST STRUCTURE/, 'summary: the title');
	like($s, qr/resolution\s+1\.85/, 'summary: the resolution');
	like($s, qr/chain A\s+protein/, 'summary: a line per chain');
	like($s, qr/MAGCMHHSC/, 'summary: the sequence');
	like($s, qr/1 gap\b/, 'summary: gaps, singular');
	like($s, qr/NAG_A_201/, 'summary: the ligands');
	is(substr($s, -1), "\n", 'summary: ends with a newline');
}

#--------
# the views all insist on a real structure
#--------
for my $f (qw(structure_atoms structure_residues structure_ligands structure_sequences structure_summary)) {
	no strict 'refs';
	# structure_sequences also takes a file name, so its complaint about being
	# handed nothing says so; every one of them still refuses a hash that is
	# not a structure.
	my $nothing = $f eq 'structure_sequences'
		? qr/expected a file name or the hash reference/
		: qr/expected the hash reference/;
	throws_ok { $f->({ not => 'a structure' }) } qr/expected the hash reference/,
		"$f: a hash that is not a structure dies";
	throws_ok { $f->(undef) } $nothing, "$f: undef dies";
	throws_ok { $f->([]) } qr/expected the hash reference/, "$f: an arrayref dies";
	throws_ok { $f->('') } $nothing, "$f: the empty string dies";
}

#--------
# h() -- the documentation, by name, by glob and by reference.
#
# h() prints to STDOUT by name rather than to the selected handle, on purpose:
# help that goes wherever output happens to have been redirected is help you
# cannot find.  So the capture has to replace the handle itself.
#--------
sub capture_h {
	my @args = @_;
	my ($out, $ret) = ('');
	{
		local *STDOUT;
		open STDOUT, '>', \$out or die $!;
		$ret = h(@args);
		close STDOUT;
	}
	return ($ret, $out);
}

{
	my ($got, $out) = capture_h('structure_info');
	is($got, 'structure_info', 'h: returns the name it showed');
	like($out, qr/structure_info/, 'h: prints that function-s documentation');
	like($out, qr/altloc/, 'h: including its options');
	unlike($out, qr/=head/, 'h: without the POD markup around it');
}
for my $arg ('res_type', \&Chem::Structure::Parser::res_type, *Chem::Structure::Parser::res_type) {
	my ($got, $out) = capture_h($arg);
	like($out, qr/nucleotide/, 'h: takes a name, a reference or a glob');
}
{
	# a fully qualified name, which is what a glob stringifies to
	my ($got, $out) = capture_h('Chem::Structure::Parser::aa3to1');
	is($got, 'aa3to1', 'h: a package-qualified name is accepted');
}
{
	my ($got, $out) = capture_h();
	is($got, undef, 'h: with no argument returns undef');
	like($out, qr/Documented functions/, 'and lists what is documented');
	like($out, qr/\bres1\b/, 'including the small ones');
}
{
	my ($got, $out) = capture_h('no_such_function');
	is($got, undef, 'h: an undocumented name returns undef');
	like($out, qr/Documented functions/, 'and lists the documented ones instead of dying');
}

done_testing();
