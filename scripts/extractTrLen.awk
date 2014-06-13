{
print $0; 
print $4; 
print gtf; 
tr_len=0;
#("grep " $4 " " gtf) | getline u | tr_len+=u;
#If you have more than one line, you can loop through the results:
#Sursa: http://www.grymoire.com/Unix/Awk.html
i=0;
print "While..";
while ( "grep " $4 " " gtf | getline ) {
		
		print $0;
		print $5;
		print $4;
		
		#cmd[i++] = $5-$4+1;
		cmd[i++] = $0;
}

for (i in cmd) {
    printf("%s=%s\n", i, cmd[i]);
}

#printf("Transcript length: %s\n", tr_len); 

}
