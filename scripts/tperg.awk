#Print number of transcripts per gene

BEGIN{
	per["none"]=0;
	#per[gene_name]
} 
{

current_gene=$10;
#remove quotes and semicolon
gsub("\"","", current_gene);
gsub(";","", current_gene);

	if( "none" in per )
	{
		#print "First Gene";
		delete per["none"];
			
		per[current_gene]=1;
		previous_tr=$12;
		previous_gene=current_gene;
	}
	else if (current_gene in per)
	{
		#print "Same gene";
		if ( $12 == previous_tr )
		{
			#print "\tSame transcript"
		}
		else
		{
			#print "\tNew Transcript"
			per[current_gene]++;
		}	
		previous_tr=$12;
		previous_gene=current_gene;
	}
	else
	{
		#print "new gene -- print the number of transcripts in the previous gene";
		print previous_gene, per[previous_gene];
		
		#Assignt the new values
		per[current_gene]=1;
		previous_gene=current_gene;
		previous_tr=$12;
	}
}
END{#Print the last gene
print previous_gene, per[previous_gene];
}