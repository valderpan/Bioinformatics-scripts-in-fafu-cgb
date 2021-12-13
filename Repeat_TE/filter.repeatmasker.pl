open IN,"$ARGV[0]" or die; #repeatmasker output
open OUT,">$ARGV[1]";
open dele,">$ARGV[1].deletion";
while (<IN>){
	next if (/^    SW  | score  /);
	chomp;$_=~s/\r//;
	if (/^\h*(.+)/){
		$item=$1;$item=~s/\h+/\t/g;my @a=split/\t/,$item;
		if ($a[4] eq $last_chr){
			if ($a[5] > $last_e){
				if ($overlap){
					$list_range="$last_s	$last_e";
					for my $temp(sort {$score{$b} <=> $score{$a}} keys %score){
						#print $temp;
						#print "$temp	$score{$temp}\n";
						(my $temp_s,my $temp_e)=$temp=~/	(\S*)	(\S*$)/;
						if ($list_range){
							#print "	$list_range";
							$list_range=~s/\s+$//;my @range=split/\t/,$list_range;undef $list_range;
							for my $order (0..$#range){
								if ($order%2 == 0){
									my $r1=$range[$order];my $r2=$range[$order+1];
									if ($temp_s>=$r1 && $temp_e<=$r2){
										$left_remaining=$temp_s-1;$right_remaining=$temp_e+1;
										if ($left_remaining-$r1+1 >= 50){$list_range.="$r1	$left_remaining	";}
										if ($r2-$right_remaining+1 >= 50){$list_range.="$right_remaining	$r2	";}
										$tran{$temp}="$temp_s	$temp_e";
									}
									elsif ($temp_s < $r1 && $temp_e>=$r1 && $temp_e <=$r2){
										$right_remaining=$temp_e+1;
										if ($r2-$right_remaining+1 >= 50){$list_range.="$right_remaining	$r2	";}
										$tran{$temp}.="$r1	$temp_e\|";
									}
									elsif ($temp_e > $r2 && $temp_s>=$r1 && $temp_s <=$r2){
										$left_remaining=$temp_s-1;
										if ($left_remaining-$r1+1 >= 50){$list_range.="$r1	$left_remaining	";}
										$tran{$temp}.="$temp_s	$r2\|";
									}
									elsif ($temp_e > $r2 && $temp_s < $r1 ){
										$tran{$temp}.="$r1	$r2\|";
									}
									else{
										$list_range.="$r1	$r2	";
									}
								}
							}
						}
						else{
							$tran{$temp}="deletion";#print "	no_list";
						}
						#$tran{$temp} ? (print "	new: $tran{$temp}\n") : (print "	new: deletion outside\n");
						$tran{$temp}="deletion" if (!$tran{$temp});
					}
					for my $temp(sort {$sort{$a} <=> $sort{$b}} keys %sort){
						my @temp=split/\h+/,$item{$temp};
						#if (!$tran{$temp}){print "$temp	without	new position\n";}
						if ($tran{$temp} eq "deletion"){
							#print OUT "$item{$temp}	deletion\n";
							print dele "$item{$temp}	deletion\n";
						}
						else{
							$tran{$temp}=~s/\|$//;
							my @temp1=split/\|/,$tran{$temp};
							if ($#temp1 == 0){
								(my $temp_s,my $temp_e)=$temp1[0]=~/(\S*)	(\S*)$/;
								$temp[5]=$temp_s;$temp[6]=$temp_e;
								print OUT join ("\t",@temp);print OUT "\n";
							}
							else{
								for my $temp11(@temp1){
									(my $temp_s,my $temp_e)=$temp11=~/(\S*)	(\S*)$/;
									$temp[5]=$temp_s;$temp[6]=$temp_e;
									print OUT join ("\t",@temp);print OUT "	multiple_split\n";
								}
							}
						}
					}
					undef %score;undef %sort;
				}
				else{
					print OUT "$last_item\n";
				}
				$last_item=$item;$last_chr=$a[4];$last_s=$a[5];$last_e=$a[6];$last_score=$a[0];
				$overlap=0;
			}
			else{
				if (!$overlap){
					$item{"$last_chr	$last_s	$last_e"}=$last_item;
					$sort{"$last_chr	$last_s	$last_e"}=$last_s;
					$score{"$last_chr	$last_s	$last_e"}=$last_score;
				}
				if ($a[6] > $last_e){$last_e=$a[6];}
				$item{"$a[4]	$a[5]	$a[6]"}=$item;
				$sort{"$a[4]	$a[5]	$a[6]"}=$a[5];
				$score{"$a[4]	$a[5]	$a[6]"}=$a[0];
				$overlap=1;
			}
		}
		else{
			if (!$last_chr){
				$last_item=$item;$last_chr=$a[4];$last_s=$a[5];$last_e=$a[6];$last_score=$a[0];
			}
			else{
				if ($overlap){
					$list_range="$last_s	$last_e";
					for my $temp(sort {$score{$b} <=> $score{$a}} keys %score){
						#print $temp;
						#print "$temp	$score{$temp}\n";
						(my $temp_s,my $temp_e)=$temp=~/	(\S*)	(\S*$)/;
						if ($list_range){
							#print "	$list_range";
							$list_range=~s/\s+$//;my @range=split/\t/,$list_range;undef $list_range;
							for my $order (0..$#range){
								if ($order%2 == 0){
									my $r1=$range[$order];my $r2=$range[$order+1];
									if ($temp_s>=$r1 && $temp_e<=$r2){
										$left_remaining=$temp_s-1;$right_remaining=$temp_e+1;
										if ($left_remaining-$r1+1 >= 50){$list_range.="$r1	$left_remaining	";}
										if ($r2-$right_remaining+1 >= 50){$list_range.="$right_remaining	$r2	";}
										$tran{$temp}="$temp_s	$temp_e";
									}
									elsif ($temp_s < $r1 && $temp_e>=$r1 && $temp_e <=$r2){
										$right_remaining=$temp_e+1;
										if ($r2-$right_remaining+1 >= 50){$list_range.="$right_remaining	$r2	";}
										$tran{$temp}.="$r1	$temp_e\|";
									}
									elsif ($temp_e > $r2 && $temp_s>=$r1 && $temp_s <=$r2){
										$left_remaining=$temp_s-1;
										if ($left_remaining-$r1+1 >= 50){$list_range.="$r1	$left_remaining	";}
										$tran{$temp}.="$temp_s	$r2\|";
									}
									elsif ($temp_e > $r2 && $temp_s < $r1 ){
										$tran{$temp}.="$r1	$r2\|";
									}
									else{
										$list_range.="$r1	$r2	";
									}
								}
							}
						}
						else{
							$tran{$temp}="deletion";#print "	no_list";
						}
						#$tran{$temp} ? (print "	new: $tran{$temp}\n") : (print "	new: deletion outside\n");
						$tran{$temp}="deletion" if (!$tran{$temp});
					}
					for my $temp(sort {$sort{$a} <=> $sort{$b}} keys %sort){
						my @temp=split/\h+/,$item{$temp};
						#if (!$tran{$temp}){print "$temp	without	new position\n";}
						if ($tran{$temp} eq "deletion"){
							#print OUT "$item{$temp}	deletion\n";
							print dele "$item{$temp}	deletion\n";
						}
						else{
							$tran{$temp}=~s/\|$//;
							my @temp1=split/\|/,$tran{$temp};
							if ($#temp1 == 0){
								(my $temp_s,my $temp_e)=$temp1[0]=~/(\S*)	(\S*)$/;
								$temp[5]=$temp_s;$temp[6]=$temp_e;
								print OUT join ("\t",@temp);print OUT "\n";
							}
							else{
								for my $temp11(@temp1){
									(my $temp_s,my $temp_e)=$temp11=~/(\S*)	(\S*)$/;
									$temp[5]=$temp_s;$temp[6]=$temp_e;
									print OUT join ("\t",@temp);print OUT "	multiple_split\n";
								}
							}
						}
					}
					undef %score;undef %sort;
				}
				else{
					print OUT "$last_item\n";
				}
				$last_item=$item;$last_chr=$a[4];$last_s=$a[5];$last_e=$a[6];$last_score=$a[0];
				$overlap=0;
			}
		}
	}
}


if ($overlap){
	$list_range="$last_s	$last_e";
	for my $temp(sort {$score{$b} <=> $score{$a}} keys %score){
		#print $temp;
		#print "$temp	$score{$temp}\n";
		(my $temp_s,my $temp_e)=$temp=~/	(\S*)	(\S*$)/;
		if ($list_range){
			#print "	$list_range";
			$list_range=~s/\s+$//;my @range=split/\t/,$list_range;undef $list_range;
			for my $order (0..$#range){
				if ($order%2 == 0){
					my $r1=$range[$order];my $r2=$range[$order+1];
					if ($temp_s>=$r1 && $temp_e<=$r2){
						$left_remaining=$temp_s-1;$right_remaining=$temp_e+1;
						if ($left_remaining-$r1+1 >= 50){$list_range.="$r1	$left_remaining	";}
						if ($r2-$right_remaining+1 >= 50){$list_range.="$right_remaining	$r2	";}
						$tran{$temp}="$temp_s	$temp_e";
					}
					elsif ($temp_s < $r1 && $temp_e>=$r1 && $temp_e <=$r2){
						$right_remaining=$temp_e+1;
						if ($r2-$right_remaining+1 >= 50){$list_range.="$right_remaining	$r2	";}
						$tran{$temp}.="$r1	$temp_e\|";
					}
					elsif ($temp_e > $r2 && $temp_s>=$r1 && $temp_s <=$r2){
						$left_remaining=$temp_s-1;
						if ($left_remaining-$r1+1 >= 50){$list_range.="$r1	$left_remaining	";}
						$tran{$temp}.="$temp_s	$r2\|";
					}
					elsif ($temp_e > $r2 && $temp_s < $r1 ){
						$tran{$temp}.="$r1	$r2\|";
					}
					else{
						$list_range.="$r1	$r2	";
					}
				}
			}
		}
		else{
			$tran{$temp}="deletion";#print "	no_list";
		}
		#$tran{$temp} ? (print "	new: $tran{$temp}\n") : (print "	new: deletion outside\n");
		$tran{$temp}="deletion" if (!$tran{$temp});
	}
	for my $temp(sort {$sort{$a} <=> $sort{$b}} keys %sort){
		my @temp=split/\h+/,$item{$temp};
		#if (!$tran{$temp}){print "$temp	without	new position\n";}
		if ($tran{$temp} eq "deletion"){
			#print OUT "$item{$temp}	deletion\n";
			print dele "$item{$temp}	deletion\n";
		}
		else{
			$tran{$temp}=~s/\|$//;
			my @temp1=split/\|/,$tran{$temp};
			if ($#temp1 == 0){
				(my $temp_s,my $temp_e)=$temp1[0]=~/(\S*)	(\S*)$/;
				$temp[5]=$temp_s;$temp[6]=$temp_e;
				print OUT join ("\t",@temp);print OUT "\n";
			}
			else{
				for my $temp11(@temp1){
					(my $temp_s,my $temp_e)=$temp11=~/(\S*)	(\S*)$/;
					$temp[5]=$temp_s;$temp[6]=$temp_e;
					print OUT join ("\t",@temp);print OUT "	multiple_split\n";
				}
			}
		}
	}
}
else{
	print OUT "$last_item\n";
}